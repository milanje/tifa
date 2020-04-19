#!/bin/sh

#
# Clean and reset the SCons environment to build the project
# from a clean state.
#
scons -c

if [ -d ./build ];
then
    cmd='rm -fr ./build'
    echo $cmd
    eval $cmd
fi

if [ -d ./.sconf_temp ];
then
    cmd='rm -fr ./.sconf_temp'
    echo $cmd
    eval $cmd
fi

if [ -f ./.sconsign.dblite ];
then
    cmd='rm -f ./.sconsign.dblite'
    echo $cmd
    eval $cmd
fi

if [ -f ./config.log ];
then
    cmd='rm -f ./config.log'
    echo $cmd
    eval $cmd
fi

if [ -f ./tifa_config.h ];
then
    cmd='rm -f ./tifa_config.h'
    echo $cmd
    eval $cmd
fi

