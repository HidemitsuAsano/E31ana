#! /bin/sh

agconf=$(printf ag:~/k18ana/conf/Run69/*)
ccconf=$(printf ./conf/Run69/)

agparam=$(printf ag:~/k18ana/param/Run69/*)
ccparam=$(printf ./param/Run69/)

rsync -a "$agconf" "$ccconf"
rsync -a "$agparam" "$ccparam"

