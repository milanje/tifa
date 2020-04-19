# Copyright (C) 2011 CNRS - Ecole Polytechnique - INRIA.
#
# This file is part of TIFA.
#
# TIFA is free software; you can redistribute it and/or modify it under the
# terms of the GNU Lesser General Public License as published by the Free
# Software Foundation; either version 2.1 of the License, or (at your option)
# any later version.
#
# TIFA is distributed in the hope that it will be useful, but WITHOUT ANY
# WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
# FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for
# more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with this library; if not, write to the Free Software Foundation,
# Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301 USA.

#------------------------------------------------------------------------------
#                        Top level SCons input file.
#------------------------------------------------------------------------------
# File    : SConstruct
# Author  : Jerome Milan
# Version : 2020-04-15
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
# Imported Python modules
#------------------------------------------------------------------------------
import os
import sys
import re

from math import sqrt
from math import log10

#------------------------------------------------------------------------------
# Command line options and default values : see BuildOptions.py for a more
# detailled description
#------------------------------------------------------------------------------
devel_conf_file = 'LocalBuildOptions.py'
devel_build     = False
if os.path.exists(os.path.join(os.getcwd(), devel_conf_file)):
    #
    # Check for a LocalBuildOptions.py file. That way, developpers can have
    # their own option file without polluting the default BuildOptions.py file
    # that will be distributed...
    #
    opts = Variables(devel_conf_file)
    devel_build = True
else:
    opts = Variables('BuildOptions.py')
#
# Installation directories options
#
installdirs = []

opts.Add(
    PathVariable(
        'PREFIX',
        'Prefix for the other intallation directories',
        '/usr/local',
        PathVariable.PathIsDirCreate
    )
)

opts.Add(
    'BINDIR',
    'Where to install programs (prefixed by PREFIX if not an absolute path)',
    'bin',
)
installdirs.append('BINDIR')

opts.Add(
    'INCLUDEDIR',
    "Where to install headers (prefixed by PREFIX if not an absolute path)",
    'include/tifa',
)
installdirs.append('INCLUDEDIR')

opts.Add(
    'SCRIPTDIR',
    'Where to install scripts (prefixed by PREFIX if not an absolute path)',
    'bin'
)
installdirs.append('SCRIPTDIR')

#
# Find standard perl lib location
#
command        = 'perl -MConfig -le \'print \"$Config{installsitelib}\"\''
std_perlmoddir = os.popen(command).read()

opts.Add(
    'PERLMODDIR',
    'Where to install perl modules (prefixed by PREFIX if not an absolute path)',
    std_perlmoddir
)
installdirs.append('PERLMODDIR')

opts.Add(
    'CONFDIR',
    'Where to install programs\' configuration files (prefixed by PREFIX if not an absolute path)',
    'share/tifa/conf'
)
installdirs.append('CONFDIR')

opts.Add(
    'LIBDIR',
    'Where to install the TIFA library (prefixed by PREFIX if not an absolute path)',
    'lib'
)
installdirs.append('LIBDIR')

opts.Add(
    'MANDIR',
    'Where to install man pages',
    'man'
)
installdirs.append('MANDIR')

opts.Add(
    'DOCDIR',
    'Where to install TIFA\'s Doxygen documentation (prefixed by PREFIX if not an absolute path)',
    'share/tifa/doc'
)
installdirs.append('DOCDIR')

#
# Doxygen documentation options
#
opts.Add(
    EnumVariable(
        'GENERATE_DOXYGEN_DOC',
        'Set to build doxygen documentation',
        'yes',
        allowed_values = ('yes', 'no'),
        ignorecase = 2
    )
)
opts.Add(
    EnumVariable(
        'DOXYGEN_PAPER_TYPE',
        'Set the paper type of the PDF Doxygen documentation',
        'a4wide',
        allowed_values=('a4', 'a4wide', 'letter', 'legal', 'executive'),
        ignorecase = 2
    )
)
opts.Add(
    EnumVariable(
        'GENERATE_DOXYGEN_HTML',
        'Set to build doxygen documentation as HTML files',
        'yes',
        allowed_values = ('yes', 'no'),
        ignorecase = 2
    )
)
opts.Add(
    EnumVariable(
        'GENERATE_DOXYGEN_PDF',
        'Set to build doxygen documentation as a PDF file',
        'yes',
        allowed_values = ('yes', 'no'),
        ignorecase = 2
    )
)
opts.Add(
    EnumVariable(
        'GENERATE_DOXYGEN_XML',
        'Set to build doxygen documentation as XML files',
        'no',
        allowed_values = ('yes', 'no'),
        ignorecase = 2
    )
)
opts.Add(
    EnumVariable(
        'GENERATE_DOXYGEN_RTF',
        'Set to build doxygen documentation as a RTF file',
        'no',
        allowed_values = ('yes', 'no'),
        ignorecase = 2
    )
)
opts.Add(
    EnumVariable(
        'GENERATE_MAN_PAGES',
        'Set to build man pages',
        'yes',
        allowed_values = ('yes', 'no'),
        ignorecase = 2
    )
)
#
# Build options
#
opts.Add(
    'CC',
    'The C compiler'
)
opts.Add(
    'CCFLAGS',
    'General options passed to the C and C++ compilers',
    '-O3'
)
opts.Add(
    'LINK',
    'The linker'
)
opts.Add(
    'LINKFLAGS',
    'General options passed to the linker'
)
opts.Add(
    'BUILDDIR',
    'Where to build the TIFA library and associated programs',
    'build'
)
opts.Add(
    BoolVariable(
        'BUILDDIR_DEPENDS_ON_PLATFORM',
        "Set to 'yes' to build in platform dependent subdirectory of BUILDDIR",
        0
    )
)
opts.Add(
    'GMP_LIBDIR',
    'Directory where the GMP library is installed',
    ''
)
opts.Add(
    'GMP_INCDIR',
    'Directory where the GMP headers are installed',
    ''
)
opts.Add(
    BoolVariable(
        'USE_GMP_INTERNAL_FUNCS',
        "Set to 'yes' to use GMP's internal functions",
        0
    )
)
opts.Add(
    'GMP_INTERNAL_INCDIR',
    'Directory where the internal GMP headers are installed',
    ''
)
opts.Add(
    BoolVariable(
        'USE_OWN_RAND',
        "Set to 'yes' to use TIFA's basic pseudo-random generator",
        0
    )
)
opts.Add(
    'BITSTRING_T',
    'Native type used to represent a string of bits',
    'uint64_t'
)
#
# Machine / architecture dependant stuff
#
opts.Add(
    EnumVariable(
        'ENDIANNESS',
        'Endianness of the processor',
        'infer',
        allowed_values = ('big', 'little', 'infer')
    )
)
opts.Add(
    'WORDSIZE',
    'Wordsize of the processor (in bits)',
    ''
)
opts.Add(
    BoolVariable(
        'USE_CALLOC_MEMSET',
        "Set to 'yes' to use calloc/memset for integer array initialization",
        1
    )
)
#
# Verbosity options
#
opts.Add(
    BoolVariable(
        'PRINT_ERROR',
        "Set to 'yes' to print critical error messages on stderr",
        1
    )
)
opts.Add(
    BoolVariable(
        'ALLOW_VERBOSE',
        "Set to 'yes' to turn on output messages in some tifa functions",
        0
    )
)
opts.Add(
    BoolVariable(
        'VERBOSE_CFRAC',
        "Set to 'yes' to turn on output messages in the cfrac function",
        1
    )
)
opts.Add(
    BoolVariable(
        'VERBOSE_ECM',
        "Set to 'yes' to turn on output messages in the ecm function",
        1
    )
)
opts.Add(
    BoolVariable(
        'VERBOSE_FERMAT',
        "Set to 'yes' to turn on output messages in the fermat function",
        1
    )
)
opts.Add(
    BoolVariable(
        'VERBOSE_SIQS',
        "Set to 'yes' to turn on output messages in the siqs function",
        1
    )
)
opts.Add(
    BoolVariable(
        'VERBOSE_SQUFOF',
        "Set to 'yes' to turn on output messages in the squfof function",
        1
    )
)
opts.Add(
    BoolVariable(
        'VERBOSE_RHO',
        "Set to 'yes' to turn on output messages in the rho function",
        1
    )
)
opts.Add(
    BoolVariable(
        'VERBOSE_TDIV',
        "Set to 'yes' to turn on output messages in the tdiv function",
        1
    )
)
#
# Timing options
#
opts.Add(
    BoolVariable(
        'ALLOW_TIMING',
        "Set to 'yes' to turn on timing measurements in some tifa functions",
        0
    )
)
opts.Add(
    BoolVariable(
        'ALLOW_EXTRA_TIMING',
        "Set to 'yes' to turn on extra timing measurements in some "
        "tifa functions",
        0
    )
)
opts.Add(
    BoolVariable(
        'TIMING_CFRAC',
        "Set to 'yes' to turn on timing messages for the cfrac function",
        1
    )
)
opts.Add(
    BoolVariable(
        'TIMING_ECM',
        "Set to 'yes' to turn on timing messages for the ecm function",
        1
    )
)
opts.Add(
    BoolVariable(
        'TIMING_FERMAT',
        "Set to 'yes' to turn on timing messages for the fermat function",
        1
    )
)
opts.Add(
    BoolVariable(
        'TIMING_SIQS',
        "Set to 'yes' to turn on timing messages for the siqs function",
        1
    )
)
opts.Add(
    BoolVariable(
        'TIMING_SQUFOF',
        "Set to 'yes' to turn on timing messages for the squfof function",
        1
    )
)
opts.Add(
    BoolVariable(
        'TIMING_RHO',
        "Set to 'yes' to turn on timing messages for the rho function",
        1
    )
)
opts.Add(
    BoolVariable(
        'TIMING_TDIV',
        "Set to 'yes' to turn on timing messages for the tdiv function",
        1
    )
)

#------------------------------------------------------------------------------
# Create a custom environement
#------------------------------------------------------------------------------
#
# Import the whole external environement
#
env = Environment(
          options = opts,
          tools   = ['default'],
          ENV     = os.environ
      )

#
# Check install and build directories
#
for directory in installdirs:
    if not os.path.isabs(directory):
        env[directory] = os.path.join(env['PREFIX'], env[directory])

if not os.path.isabs(env['BUILDDIR']):
    env['BUILDDIR'] = os.path.join(os.getcwd(), env['BUILDDIR'])

if env['BUILDDIR_DEPENDS_ON_PLATFORM'] != 0:
    print("Infering build directory... ")
    try:
        import platform
        system = platform.system().lower()
        proc   = platform.processor().lower()
        env['BUILDDIR'] = os.path.join(env['BUILDDIR'], proc, system)
    except ImportError:
        print("Warning: Cannot find python module 'platform':")
        print("         -> platform build directory based on OS name only.")
        platform = sys.platform.lower()
        env['BUILDDIR'] = os.path.join(env['BUILDDIR'], platform)
    except:
        print("Warning: Cannot infer build directory:")
        print("         -> using directory BUILDDIR from BuildOptions.py.")

sysinfo = os.uname();

env['COMPILED_ON_OS']      = sysinfo[0];
env['COMPILED_ON_RELEASE'] = sysinfo[2];
env['COMPILED_ON_MACHINE'] = sysinfo[4];
env['COMPILED_ON']         = '-'.join([sysinfo[0], sysinfo[2],
                                       sysinfo[4], env['WORDSIZE'] + "bit"])

#
# Check verbosity options
#
if env['ALLOW_VERBOSE'] == 0:
    env['VERBOSE_CFRAC']  = 0
    env['VERBOSE_ECM']    = 0
    env['VERBOSE_FERMAT'] = 0
    env['VERBOSE_SIQS']   = 0
    env['VERBOSE_SQUFOF'] = 0
    env['VERBOSE_RHO']    = 0
    env['VERBOSE_TDIV']   = 0
#
# Check timing measurements options
#
if env['ALLOW_TIMING'] == 0:
    env['ALLOW_EXTRA_TIMING'] = 0
    env['TIMING_CFRAC']       = 0
    env['TIMING_ECM']         = 0
    env['TIMING_FERMAT']      = 0
    env['TIMING_SIQS']        = 0
    env['TIMING_SQUFOF']      = 0
    env['TIMING_RHO']         = 0
    env['TIMING_TDIV']        = 0
#
# GMP related variables
#
if (env['GMP_LIBDIR'] != ""):
    gmplibdir = Dir(env['GMP_LIBDIR'])
    env.Prepend(LIBPATH = [str(gmplibdir)])

if (env['GMP_INCDIR'] != ""):
    gmpincdir = Dir(env['GMP_INCDIR'])
    env.Prepend(CPPPATH = [str(gmpincdir)])

if (env['USE_GMP_INTERNAL_FUNCS'] != 0):
    if (env['GMP_INTERNAL_INCDIR'] != ""):
        gmpincdir = Dir(env['GMP_INTERNAL_INCDIR']).abspath
        env.Prepend(CPPPATH = [str(gmpincdir)])

#------------------------------------------------------------------------------
# Generate help text given by "scons -h"
#------------------------------------------------------------------------------
Help(opts.GenerateHelpText(env))

#------------------------------------------------------------------------------
# Configure for host: checks for libraries, header files, functions, etc...
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
# Some checking subroutines
#------------------------------------------------------------------------------
def checkSizeOf(context, type):
    context.Message( 'Getting size of ' + type + '... ' )
    #
    # We have to be careful with the program we wish to test here since
    # compilation will be attempted using the current environment's flags.
    # So make sure that the program will compile without any warning. For
    # example using: 'int main(int argc, char** argv)' will fail with the
    # '-Wall -Werror' flags since the variables argc and argv would not be
    # used in the program...
    #
    program = """
#include <stdlib.h>
#include <stdint.h>
#include <stdbool.h>
#include <stdio.h>
int main() {
    printf("%d", (int)sizeof(""" + type + """));
    return 0;
}
"""
    ret = context.TryCompile(program, '.c')
    if ret == 0:
        print("test failed!\n")
        print("Error: cannot run the following test program:")
        print(program)
        print("Please check your compiler flags.")
        Exit(1)
    ret = context.TryRun(program, '.c')
    context.Result(ret[1])
    return ret[1]

#------------------------------------------------------------------------------
def checkSymbolValue(context, symbol, fmt):
    context.Message( 'Getting value of symbol ' + symbol + '... ' )
    program = """
#include <stdlib.h>
#include <stdio.h>
#include <gmp.h>

int main() {
  printf(\"""" + fmt + """\", """ + symbol + """);
  return 0;
}
"""
    ret = context.TryCompile(program, '.c')
    if ret == 0:
        print("test failed!\n")
        print("Error: cannot run the following test program:")
        print(program)
        print("Please check your compiler flags.")
        Exit(1)
    ret = context.TryRun(program, '.c')
    context.Result(ret[1])
    return ret[1]

#------------------------------------------------------------------------------
def checkCharBit(context):
    context.Message( 'Getting number of bits in a char... ' )
    program = """
#include <stdlib.h>
#include <stdio.h>
#include <limits.h>
int main() {
    printf("%d", (int)CHAR_BIT);
    return 0;
}
"""
    ret = context.TryCompile(program, '.c')
    if ret == 0:
        print("test failed!\n")
        print("Error: cannot run the following test program:")
        print(program)
        print("Please check your compiler flags.")
        Exit(1)
    ret = context.TryRun(program, '.c')
    context.Result(ret[1])
    return ret[1]

#------------------------------------------------------------------------------
def checkUlongMax(context):
    context.Message( 'Getting ULONG_MAX... ' )
    program = """
#include <stdlib.h>
#include <stdio.h>
#include <limits.h>
int main() {
    printf("%lu", (unsigned long int)ULONG_MAX);
    return 0;
}
"""
    ret = context.TryCompile(program, '.c')
    if ret == 0:
        print("test failed!\n")
        print("Error: cannot run the following test program:")
        print(program)
        print("Please check your compiler flags.")
        Exit(1)
    ret = context.TryRun(program, '.c')
    context.Result(ret[1])
    return ret[1]

#------------------------------------------------------------------------------
def checkLongMax(context):
    context.Message( 'Getting LONG_MAX... ' )
    program = """
#include <stdlib.h>
#include <stdio.h>
#include <limits.h>
int main() {
    printf("%lu", (unsigned long int)LONG_MAX);
    return 0;
}
"""
    ret = context.TryCompile(program, '.c')
    if ret == 0:
        print("test failed!\n")
        print("Error: cannot run the following test program:")
        print(program)
        print("Please check your compiler flags.")
        Exit(1)
    ret = context.TryRun(program, '.c')
    context.Result(ret[1])
    return ret[1]
#------------------------------------------------------------------------------
def checkEndianness(context):
    context.Message("Checking for endianness... ")
    if context.env['ENDIANNESS'] == "infer":
        try:
            context.Result(sys.byteorder)
            return str(sys.byteorder)
        except:
            print("ERROR: Cannot infer endianness:")
            print("       -> Set ENDIANNESS manually in BuildOptions.py.")
            sys.exit(-1)
#------------------------------------------------------------------------------
def checkSizeWord(context):
    if (context.env['WORDSIZE'] == "") or (context.env['WORDSIZE'] == "infer"):
        context.Message("Checking for size of machine word... ")
        try:
            import platform
            res = platform.architecture()[0][0:2]
            context.Result(res)
            return res
        except ImportError:
            print("Warning: Cannot find python module 'platform':")
            print("         -> wordsize based on max integer.")
            if sys.maxint == 65535:
                context.Result(16)
                return 16
            elif sys.maxint == 2147483647:
                context.Result(32)
                return 32
            elif sys.maxint == 9223372036854775807:
                context.Result(64)
                return 64
            else:
                print("Warning: Cannot infer size of machine word:")
                print("         -> Set WORDSIZE manually in BuildOptions.py.")
        except:
            print(sys.maxint)
            print("ERROR: Cannot infer size of machine word:")
            print("       -> Set WORDSIZE manually in BuildOptions.py.")
            sys.exit(-1)

#------------------------------------------------------------------------------
conf = Configure(
           env,
           custom_tests = {
               'checkSizeOf'      : checkSizeOf,
               'checkSymbolValue' : checkSymbolValue,
               'checkCharBit'     : checkCharBit,
               'checkUlongMax'    : checkUlongMax,
               'checkLongMax'     : checkLongMax,
               'checkEndianness'  : checkEndianness,
               'checkSizeWord'    : checkSizeWord
           }
       )
#
# Check for the Math library and header...
#
if not conf.CheckLib('m'):
    print('Did not find the math library. Exiting...')
    Exit(1)
if not conf.CheckHeader('math.h'):
    print('Did not find the math.h header. Exiting...')
    Exit(1)

#
# Check for the GMP library and header...
#
if not conf.CheckLib('gmp'):
    print('Did not find the GMP library. Exiting...')
    Exit(1)
if not conf.CheckHeader('gmp.h'):
    print('Did not find the gmp.h header. Exiting...')
    Exit(1)
if (env['USE_GMP_INTERNAL_FUNCS'] != 0):
    #
    # Make sure that we include gmp.h before testing for gmp-impl.h otherwise
    # SCons will not recognise gmp-impl.h as an useable header... even if the
    # file exists...
    #
    if not conf.CheckHeader(['gmp.h', 'gmp-impl.h']):
        print('Did not find the gmp-impl header. Exiting...')
        Exit(1)

#
# Check for the native type used to represent strings of bits and extract
# size for use in macros
#
conf.env['SIZEOF_BITSTRING_T'] = conf.checkSizeOf(env['BITSTRING_T'])
conf.env['CHAR_BIT']           = conf.checkCharBit()
conf.env['ULONG_MAX']          = conf.checkUlongMax()
conf.env['SIZEOF_ULONG_T']     = conf.checkSizeOf('unsigned long int')
conf.env['LONG_MAX']           = conf.checkLongMax()
conf.env['SIZEOF_LONG_T']      = conf.checkSizeOf('long int')

#
# Machine or architecture dependant checks
#
conf.env['ENDIANNESS']   = conf.checkEndianness()
conf.env['WORDSIZE']     = conf.checkSizeWord()
conf.env['GMP_NUMB_MAX'] = conf.checkSymbolValue('GMP_NUMB_MAX', '%lu')

env = conf.Finish()

env['BITSTRING_T_BITSIZE'] = int(env['SIZEOF_BITSTRING_T']) \
                               * int(env['CHAR_BIT'])
env['ULONG_T_BITSIZE']     = int(env['SIZEOF_ULONG_T']) \
                               * int(env['CHAR_BIT'])
env['ULONG_MAX_DIVIDED_BY_2'] = int(env['ULONG_MAX']) // 2
env['ULONG_MAX_DIVIDED_BY_3'] = int(env['ULONG_MAX']) // 3
env['ULONG_MAX_DIVIDED_BY_4'] = int(env['ULONG_MAX']) // 4
env['ULONG_MAX_DIVIDED_BY_5'] = int(env['ULONG_MAX']) // 5
env['ULONG_MAX_DIVIDED_BY_6'] = int(env['ULONG_MAX']) // 6
env['SQRT_ULONG_MAX']         = int(sqrt(int(env['ULONG_MAX'])))

env['LONG_T_BITSIZE']        = int(env['SIZEOF_LONG_T']) \
                                 * int(env['CHAR_BIT'])
env['LONG_MAX_DIVIDED_BY_2'] = int(env['LONG_MAX']) // 2
env['LONG_MAX_DIVIDED_BY_3'] = int(env['LONG_MAX']) // 3
env['LONG_MAX_DIVIDED_BY_4'] = int(env['LONG_MAX']) // 4
env['LONG_MAX_DIVIDED_BY_5'] = int(env['LONG_MAX']) // 5
env['LONG_MAX_DIVIDED_BY_6'] = int(env['LONG_MAX']) // 6
env['SQRT_LONG_MAX']         = int(sqrt(int(env['LONG_MAX'])))
env['LOG10_GMP_NUMB_MAX']    = log10(int(env['GMP_NUMB_MAX']))

#------------------------------------------------------------------------------
# Directories definitions
#------------------------------------------------------------------------------
top_dir   = os.getcwd()
lib_dir   = os.path.join(top_dir, 'lib')
tools_dir = os.path.join(top_dir, 'tools')

lib_scons_file   = os.path.join(lib_dir, 'SConscript')
tools_scons_file = os.path.join(tools_dir, 'SConscript')

base_build_dir  = os.path.join(top_dir, env['BUILDDIR'])
lib_build_dir   = os.path.join(base_build_dir, 'lib')
tools_build_dir = os.path.join(base_build_dir, 'tools')

#------------------------------------------------------------------------------
# Misc. variables
#------------------------------------------------------------------------------
env['VERSION_MAJOR']  = '0'
env['VERSION_MINOR']  = '1'
env['VERSION_PATCH']  = '4'
if devel_build:
    import time
    today = time.localtime()
    env['VERSION_DEVTAG'] = '%.4s%.2d%.2d' % (today[0], today[1], today[2])
else:
    #
    # Quick and _really_ dirty hack to tag the released development version
    #
    devel_build = 1
    env['VERSION_DEVTAG'] = '20150121'

env['VERSION'] = '.'.join([env['VERSION_MAJOR'],
                           env['VERSION_MINOR'],
                           env['VERSION_PATCH']
                          ])

if devel_build:
    env['VERSION'] += ' (devel:' + env['VERSION_DEVTAG'] + ")"

env['SHORTNAME']   = 'TIFA'
env['FULLNAME']    = "TIFA - Tools for Integer FActorization"
env['TIFA_STRING'] = env['SHORTNAME'] + ' v' + env['VERSION']
env['TARNAME']     = "tifa"

tifa_libname  = 'tifa'

tifa_distname = env['TARNAME'] + '-' + env['VERSION']
tifa_distname = re.sub("[^a-zA-Z0-9_\.]+", "_", tifa_distname);
tifa_distname = re.sub("_$", "", tifa_distname);
tifa_distname = re.sub("^_", "", tifa_distname);
tifa_distname = os.path.join(os.getcwd(), tifa_distname)

#------------------------------------------------------------------------------
# Create the project's tifa_config.h file
#------------------------------------------------------------------------------
subst_dict = {}

subst_dict['%SHORTNAME%'] = env['SHORTNAME']
subst_dict['%FULLNAME%']  = env['FULLNAME']
subst_dict['%STRING%']    = env['TIFA_STRING']
subst_dict['%TARNAME%']   = env['TARNAME']

subst_dict['%BUGREPORT%'] = 'milanje at gmail dot com'

subst_dict['%VERSION_MAJOR%'] = env['VERSION_MAJOR']
subst_dict['%VERSION_MINOR%'] = env['VERSION_MINOR']
subst_dict['%VERSION_PATCH%'] = env['VERSION_PATCH']
subst_dict['%VERSION%']       = env['VERSION']

subst_dict['%COMPILED_ON_OS%']      = env['COMPILED_ON_OS']
subst_dict['%COMPILED_ON_RELEASE%'] = env['COMPILED_ON_RELEASE']
subst_dict['%COMPILED_ON_MACHINE%'] = env['COMPILED_ON_MACHINE']
subst_dict['%COMPILED_ON%']         = env['COMPILED_ON']

subst_dict['%ENDIANNESS%'] = env['ENDIANNESS']
subst_dict['%WORDSIZE%']   = env['WORDSIZE']

subst_dict['%USE_GMP_INTERNAL_FUNCS%'] = str(int(env['USE_GMP_INTERNAL_FUNCS']))
subst_dict['%USE_OWN_RAND%'] = str(int(env['USE_OWN_RAND']))

subst_dict['%USE_CALLOC_MEMSET%']   = str(int(env['USE_CALLOC_MEMSET']))
subst_dict['%BITSTRING_T%']         = env['BITSTRING_T']
subst_dict['%SIZEOF_BITSTRING_T%']  = str(env['SIZEOF_BITSTRING_T'])
subst_dict['%BITSTRING_T_BITSIZE%'] = str(env['BITSTRING_T_BITSIZE'])

subst_dict['%SIZEOF_ULONG_T%']         = str(env['SIZEOF_ULONG_T'])
subst_dict['%ULONG_T_BITSIZE%']        = str(env['ULONG_T_BITSIZE'])
subst_dict['%ULONG_MAX%']              = env['ULONG_MAX'] + "UL"
subst_dict['%ULONG_MAX_DIVIDED_BY_2%'] = str(env['ULONG_MAX_DIVIDED_BY_2'])+"UL"
subst_dict['%ULONG_MAX_DIVIDED_BY_3%'] = str(env['ULONG_MAX_DIVIDED_BY_3'])+"UL"
subst_dict['%ULONG_MAX_DIVIDED_BY_4%'] = str(env['ULONG_MAX_DIVIDED_BY_4'])+"UL"
subst_dict['%ULONG_MAX_DIVIDED_BY_5%'] = str(env['ULONG_MAX_DIVIDED_BY_5'])+"UL"
subst_dict['%ULONG_MAX_DIVIDED_BY_6%'] = str(env['ULONG_MAX_DIVIDED_BY_6'])+"UL"
subst_dict['%SQRT_ULONG_MAX%']         = str(env['SQRT_ULONG_MAX'])+"UL"

subst_dict['%SIZEOF_LONG_T%']         = str(env['SIZEOF_LONG_T'])
subst_dict['%LONG_T_BITSIZE%']        = str(env['LONG_T_BITSIZE'])
subst_dict['%LONG_MAX%']              = env['LONG_MAX'] + "L"
subst_dict['%LONG_MAX_DIVIDED_BY_2%'] = str(env['LONG_MAX_DIVIDED_BY_2']) + "L"
subst_dict['%LONG_MAX_DIVIDED_BY_3%'] = str(env['LONG_MAX_DIVIDED_BY_3']) + "L"
subst_dict['%LONG_MAX_DIVIDED_BY_4%'] = str(env['LONG_MAX_DIVIDED_BY_4']) + "L"
subst_dict['%LONG_MAX_DIVIDED_BY_5%'] = str(env['LONG_MAX_DIVIDED_BY_5']) + "L"
subst_dict['%LONG_MAX_DIVIDED_BY_6%'] = str(env['LONG_MAX_DIVIDED_BY_6']) + "L"
subst_dict['%SQRT_LONG_MAX%']         = str(env['SQRT_LONG_MAX']) + "L"
subst_dict['%LOG10_GMP_NUMB_MAX%']    = str(env['LOG10_GMP_NUMB_MAX'])

subst_dict['%PRINT_ERROR%']    = str(int(env['PRINT_ERROR']))
subst_dict['%VERBOSE_CFRAC%']  = str(int(env['VERBOSE_CFRAC']))
subst_dict['%VERBOSE_ECM%']    = str(int(env['VERBOSE_ECM']))
subst_dict['%VERBOSE_FERMAT%'] = str(int(env['VERBOSE_FERMAT']))
subst_dict['%VERBOSE_SIQS%']   = str(int(env['VERBOSE_SIQS']))
subst_dict['%VERBOSE_SQUFOF%'] = str(int(env['VERBOSE_SQUFOF']))
subst_dict['%VERBOSE_RHO%']    = str(int(env['VERBOSE_RHO']))
subst_dict['%VERBOSE_TDIV%']   = str(int(env['VERBOSE_TDIV']))

subst_dict['%ALLOW_EXTRA_TIMING%'] = str(int(env['ALLOW_EXTRA_TIMING']))

subst_dict['%TIMING_CFRAC%']  = str(int(env['TIMING_CFRAC']))
subst_dict['%TIMING_ECM%']    = str(int(env['TIMING_ECM']))
subst_dict['%TIMING_FERMAT%'] = str(int(env['TIMING_FERMAT']))
subst_dict['%TIMING_SIQS%']   = str(int(env['TIMING_SIQS']))
subst_dict['%TIMING_SQUFOF%'] = str(int(env['TIMING_SQUFOF']))
subst_dict['%TIMING_RHO%']    = str(int(env['TIMING_RHO']))
subst_dict['%TIMING_TDIV%']   = str(int(env['TIMING_TDIV']))

config_h_in = 'tifa_config.h.in'
config_h    = 'tifa_config.h'

env.Substfile(config_h_in, SUBST_DICT = subst_dict)

#------------------------------------------------------------------------------
# Files to be distributed in the TIFA tar package
#------------------------------------------------------------------------------
env['TARFLAGS']  = '-czvf'
env['TARCOM']    = '$TAR $TARFLAGS $TARGET $SOURCES'
env['TARSUFFIX'] = '.tgz'
#
# The suffix is not automatically appended if the last dot in the filename base
# is not followed by digits only. This will be the case for the "devel" archive,
# so we make sure to append the suffix before proceeding with the Tar builder...
#
tifa_distname += env['TARSUFFIX']

if "dist" in COMMAND_LINE_TARGETS:
    archive = env.Tar(tifa_distname, 'readme')
    env.Tar(tifa_distname, 'archive.pl')
    env.Tar(tifa_distname, 'scons_files')
    env.Tar(tifa_distname, 'BuildOptions.py')
    env.Tar(tifa_distname, 'SConstruct')
    env.Tar(tifa_distname, 'doxygen.in')
    env.Tar(tifa_distname, 'tifa_config.h.in')
    env.Tar(tifa_distname, 'reset.sh')
    env.Alias('dist', archive)

#------------------------------------------------------------------------------
# Exported variables
#------------------------------------------------------------------------------
Export('env')
Export('top_dir')
Export('lib_dir')
Export('lib_build_dir')
Export('tools_dir')
Export('tifa_libname')
Export('tifa_distname')

#------------------------------------------------------------------------------
# Recurse in the lib directory to build the TIFA library
#------------------------------------------------------------------------------
SConscript(lib_scons_file, variant_dir = lib_build_dir, duplicate = 0)

#------------------------------------------------------------------------------
# Recurse in the tools directory to build the programs
#------------------------------------------------------------------------------
SConscript(tools_scons_file, variant_dir = tools_build_dir, duplicate = 0)
