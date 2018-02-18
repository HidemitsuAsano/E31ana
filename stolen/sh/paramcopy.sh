#! /bin/sh

agconf=$(printf ag:~/k18ana/conf/Run65/*)
ccconf=$(printf ./conf/Run65/)

agparam=$(printf ag:~/k18ana/param/Run65/*)
ccparam=$(printf ./param/Run65/)

rsync -a "$agconf" "$ccconf"
rsync -a "$agparam" "$ccparam"

