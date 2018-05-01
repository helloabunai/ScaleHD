from setuptools import setup, find_packages
from codecs import open
from os import path

# Get the long description from the relevant file
here = path.abspath(path.dirname(__file__))
with open(path.join(here, 'README.md'), encoding='utf-8') as f:
	long_description = f.read()

##
## Begin the setup installation
setup(
    name='ScaleHD',

    # https://packaging.python.org/en/latest/single_source_version.html
    version='0.312',

    description='Automated DNA micro-satellite genotyping.',
    long_description=long_description,
	python_requires='>=2.7.13',
    # The project's main homepage.
    url='https://github.com/helloabunai/ScaleHD',

    # Author details
    author='Alastair Maxwell/University of Glasgow',
    author_email='alastair.maxwell@glasgow.ac.uk',

    # License to ship the package with
    license='GPLv3',

    # More information?
    # See https://pypi.python.org/pypi?%3Aaction=list_classifiers
	# 2 - Pre-Alpha
	# 3 - Alpha
	# 4 - Beta
	# 5 - Stable
    classifiers=[
        'Development Status :: 2 - Pre-Alpha',

        # Indicate who your project is intended for
        'Intended Audience :: Education',
		'Intended Audience :: End Users/Desktop',
		'Intended Audience :: Science/Research',

        # Classifier matching license flag from above
        'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',

		# Specific version of the python interpreter that are supported
		# by this package. Python 3 not support at this time.
        'Programming Language :: Python :: 2.6',
        'Programming Language :: Python :: 2.7',

		## And so on
        'Environment :: Console',
		'Operating System :: MacOS :: MacOS X',
		'Operating System :: POSIX'
    ],

    # What does the project relate to?
    keywords='machine-learning genotyping bioinformatics huntington-disease',

    # You can just specify the packages manually here if your project is
    # simple. Or you can use find_packages().
    packages=find_packages(exclude=['input',
									'lib',
									'ScaleHD.egg-info',
									'build',
									'dist',
									'logs'
									]),

    # List run-time dependencies here.  These will be installed by pip when
    # your project is installed. For an analysis of "install_requires" vs pip's
    # requirements files see:
    # https://packaging.python.org/en/latest/requirements.html
    install_requires=['lxml',
					  'numpy',
					  'seaborn',
					  'matplotlib',
					  'sklearn',
					  'scipy',
					  'peakutils',
                      'pandas',
					  'pysam',
					  'regex',
					  'PyPDF2',
					  'reportlab',
					  'generatr',
					  'batchadapt',
					  'pyvcf'
					  ],

    # These are the data files to be included in the package
	# For GenoCall, this will be the data-sets used for machine-learning
	# training, and generating predictive models for each 'data state'
    package_data={'ScaleHD': ['train/long_descr.rst',
                              'train/polyglutamine.csv',
                              'train/likelihood_matrix.csv',
                              'train/raw_matrix.csv',
							  'config/config.dtd']},
	include_package_data=True,

	# Executable scripts require an entry point to allow cython to generate
	# executables for the respective target platform. This entry point is akin
	# to launching the script in bash: if __name__ == '__main__' etc..
    entry_points={
        'console_scripts': ['scalehd=ScaleHD.sherpa:main',],
    },
)