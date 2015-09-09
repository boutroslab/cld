use 5.006001;
use ExtUtils::MakeMaker;

$CC = 'c++';

# See lib/ExtUtils/MakeMaker.pm for details of how to influence
# the contents of the Makefile that is written.
WriteMakefile(
    NAME              => 'Set::IntervalTree',
    VERSION_FROM      => 'lib/Set/IntervalTree.pm', # finds $VERSION
    PREREQ_PM         => {}, # e.g., Module::Name => 1.1
    ($] >= 5.005 ?     ## Add these new keywords supported since 5.005
      (ABSTRACT_FROM  => 'lib/Set/IntervalTree.pm', # retrieve abstract from module
       AUTHOR         => 'Ben Booth <benbooth@gmail.com> and Oliver Dreier') : ()),
    LIBS              => [''], # e.g., '-lm'
    DEFINE            => '-Isrc', # e.g., '-DHAVE_SOMETHING'
    CC                => $CC,
    LD                => '$(CC)',
    OPTIMIZE          => '-g -O2',
    XSOPT             => '-C++',
    TYPEMAPS          => ['perlobject.map'],
    INC               => '-Isrc', # e.g., '-I. -I/usr/include/other'
	# Un-comment this if you add C files to link with later:
    # OBJECT            => '$(O_FILES)', # link all the C files too
);

