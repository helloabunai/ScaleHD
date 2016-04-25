from setuptools import setup, find_packages
from codecs import open
from os import path

here = path.abspath(path.dirname(__file__))

# Get the long description from the relevant file
with open(path.join(here, '..', 'README.md'), encoding='utf-8') as f:
    long_description = f.read()

setup(
    name='ScaleHD',

    # https://packaging.python.org/en/latest/single_source_version.html
    version='0.0.1',

    description='Automated DNA micro-satellite genotyping.',
    long_description=long_description,

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
        'Intended Audience :: Genetic Researchers',
		'Intended Audience :: End User/Desktop',
		'Intended Audience :: Academics',
        'Topic :: Data Analysis :: Data Classification',
		'Topic :: Machine Learning :: Support Vectors'

        # Classifier matching license flag from above
        'License :: GPLv3 :: GPLv3 License',

		# Specific version of the python interpreter that are supported
		# by this package. Python 3 not support at this time.
        'Programming Language :: Python :: 2.6',
        'Programming Language :: Python :: 2.7',

		## And so on
        'Environment :: Console',
		'Operating System :: MacOS :: MacOS X 10.8'
		'Operating System :: MacOS :: MacOS X 10.9'
		'Operating System :: MacOS :: MacOS X 10.10'
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
    install_requires=['peppercorn',
					  'progress',
					  'lxml',
					  'numpy',
					  'prettyplotlib',
					  'sklearn',
					  'scipy',
					  'pandas',
					  ],

    # These are the data files to be included in the package
	# For GenoCall, this will be the data-sets used for machine-learning
	# training, and generating predictive models for each 'data state'
    package_data={'ScaleHD': ['train/long_descr.csv', 'train/polyglu_zygosity.rst']},

	include_package_data=True,

	# Executable scripts require an entry point to allow cython to generate
	# executables for the respective target platform. This entry point is akin
	# to launching the script in bash: if __name__ == '__main__' etc..
    entry_points={
        'console_scripts': ['scalehd=ScaleHD.__sherpa:main',],
    },
)