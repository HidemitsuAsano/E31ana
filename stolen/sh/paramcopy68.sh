#! /bin/sh

agconf=$(printf ag:~/k18ana/conf/Run68/*)
ccconf=$(printf ./conf/Run68/)

agparam=$(printf ag:~/k18ana/param/Run68/*)
ccparam=$(printf ./param/Run68/)

rsync -a "$agconf" "$ccconf"
rsync -a "$agparam" "$ccparam"

