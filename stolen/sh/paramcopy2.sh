#! /bin/sh

agconf=$(printf ./conf/Run65/*)
ccconf=$(printf ~/public/ymg/k18ana/conf/Run65/)

agparam=$(printf ./param/Run65/*)
ccparam=$(printf ~/public/ymg/k18ana/param/Run65/)

cp "$agconf" "$ccconf"
cp "$agparam" "$ccparam"

