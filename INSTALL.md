# Installation #

* Get the code:
`git clone https://github.com/viyh/diphot.git`

* Install the latest numpy library locally:
`/usr/local/bin/pip2.7 install --user numpy --upgrade`

* Install the latest matplotlib library locally:
`/usr/local/bin/pip2.7 install --user matplotlib --upgrade`

* Add the local libraries to your local python search path. Edit your ~/.bash_profile and add the following:
`export PYTHONPATH=~/.local/lib/python2.7/site-packages:$PYTHONPATH`

If you are using a remote server, you will need to use SSH forwarding (`-Y`) when SSHing into the server in order to display images.
