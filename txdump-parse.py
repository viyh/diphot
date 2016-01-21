#!/usr/bin/python

#
# DiPhot: Txdump Parse - 2016-01-20
# https://github.com/viyh/diphot
#
# Parse a dump file and attempt to align stars by ID
#

# The following columns need to be dumped:
#    IMAGE,ID,XCENTER,YCENTER,OTIME,MAG,MERR

from collections import defaultdict, OrderedDict
import csv
import sys

px_threshold = 75.0
mag_threshold = 0.8
skip_px_threshold = 90.0
skip_mag_threshold = 1.0
assume = False

def read_dump(dump_filename):
    dump = defaultdict(list)
    with open(dump_filename) as dumpfile:
        for line in dumpfile.readlines():
            image, id, x, y, time, mag, merr = line.split()
            dump[image].append({'time': time, 'x': float(x), 'y': float(y), 'mag': mag, 'merr': merr})
    return OrderedDict(sorted(dump.items()))

def px_test(c1, c2):
    return abs(c1 - c2) < px_threshold

def skip_px_test(c1, c2):
    return abs(c1 - c2) > skip_px_threshold

def mag_test(star1, star2):
    if not star1.has_key('mag') or not star2.has_key('mag'):
        return True
    return star1['mag'].strip() == 'INDEF' or \
        star2['mag'].strip() == 'INDEF' or \
        abs(float(star1['mag']) - float(star2['mag'])) < mag_threshold

def skip_mag_test(star1, star2):
    if not star1.has_key('mag') or not star2.has_key('mag'):
        return False
    return star1['mag'].strip() != 'INDEF' and \
        star2['mag'].strip() != 'INDEF' and \
        abs(float(star1['mag']) - float(star2['mag'])) > skip_mag_threshold

def find_similar_star(sorted_dump, star, image):
    if not sorted_dump:
        return False
    for star_id, star_data in sorted_dump.iteritems():
        last_image, last_data = last(star_data)
        if px_test(last_data['x'], star['x']) and \
        px_test(last_data['y'], star['y']) and \
        mag_test(last_data, star):
            return star_id

    # no match, lets try manual
    print '\nimage: {}, time: {}, x: {}, y: {}, mag: {}, merr: {}'.format(image, star['time'], star['x'], star['y'], star['mag'], star['merr'])
    for star_id, star_data in sorted_dump.iteritems():
        last_image, last_data = last(star_data)
        if skip_px_test(last_data['x'], star['x']) or \
        skip_px_test(last_data['y'], star['y']) or \
        skip_mag_test(last_data, star):
            continue
        print '\timage: {}, time: {}, x: {}, y: {}, mag: {}, merr: {}'.format(last_image, last_data['time'], last_data['x'], last_data['y'], last_data['mag'], last_data['merr']),
        print "\t\tIs this the above star? ",
        if assume == False:
            continue
        if assume == True or raw_input().lower() == 'y':
            return star_id
        print '\n'
    return False

def get_new_star_id(sorted_dump):
    if sorted_dump:
        return max(sorted_dump.keys()) + 1
    else:
        return 1

def last(ordered_dict):
    key = next(reversed(ordered_dict))
    return (key, ordered_dict[key])

def sort_dump(dump):
    sorted_dump = defaultdict(OrderedDict)
    for image, data in dump.iteritems():
        for star in data:
            star_id = find_similar_star(sorted_dump, star, image)
            if not star_id:
                star_id = get_new_star_id(sorted_dump)
                print 'New star found - id: {}, image: {}, time: {}, x: {}, y: {}, mag: {}, merr: {}'.format(star_id, image, star['time'], star['x'], star['y'], star['mag'], star['merr'])
            sorted_dump[star_id][image] = star
    return sorted_dump

def write_csv(images, sorted_dump):
    with open('data.csv', 'w') as csvfile:
        csvhandle = csv.writer(csvfile, delimiter=',')
        csvhandle.writerow(['id', 'image', 'time', 'x', 'y', 'mag', 'merr'])
        for star, datapoints in sorted_dump.iteritems():
            for image in sorted(images):
                if datapoints.has_key(image):
                    datapoint = datapoints[image]
                else:
                    datapoint = {'time': 0, 'x': 0, 'y': 0, 'mag': 'INDEF', 'merr': 'INDEF'}
                row = [star, image, datapoint['time'], '{:.3f}'.format(datapoint['x']), '{:.3f}'.format(datapoint['y']), '{}'.format(datapoint['mag']), '{}'.format(datapoint['merr'])]
                csvhandle.writerow(row)

if __name__ == '__main__':
    dump_file = 'txdump.txt'
    if len(sys.argv) > 0:
        dump_file = sys.argv[1]

    dump = read_dump(dump_file)
    dump = sort_dump(dump)

    print
    images = set()
    for star in dump:
        images |= set(dump[star].keys())

    final_dump = defaultdict(OrderedDict)

    for star, values in dump.iteritems():
        star_images = set(values.keys())
        if set(images).difference(star_images):
            print "Not a complete data set: {} [missing {} datapoints out of {}]".format(star, len(set(images).difference(star_images)), len(images))
        else:
            print "Complete data set for star [{}]!".format(star)

        if float( len(set(images).difference(star_images)) ) / float ( len(images) ) * 100.0 < 3:
            print "Saving star [{}]!".format(star)
            final_dump[star] = values

    write_csv(images, final_dump)


