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
    version='1.3.0',

    description='Automated DNA micro-satellite genotyping.',
    long_description=long_description,
	long_description_content_type='text/markdown',
	python_requires='>=3.8',
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
		# by this package.
        'Programming Language :: Python :: 3.8',

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
					  'numpy==2.2.1',
					  'seaborn',
					  'matplotlib',
					  'scikit-learn==1.6.1',
					  'scipy==1.15.1',
					  'peakutils',
					  'pandas',
					  'pysam==0.22.1',
					  'regex',
					  'PyPDF2',
					  'reportlab',
					  'generatr==1.0',
					  'batchadapt==1.2',  # 1.0 was yanked
					  'cutadapt==4.9',
					  'PyVCF3',
					  'fadapa',
					  'pytest',
					  ],

    # These are the data files to be included in the package
	# For GenoCall, this will be the data-sets used for machine-learning
	# training, and generating predictive models for each 'data state'
    package_data={'ScaleHD': ['train/long_descr.rst',
                              'train/polyglutamine.csv',
                              'train/likelihood_matrix.csv',
                              'train/raw_matrix.csv',
							  'config/config.dtd',
							  'genHTML/chartBase.js',
							  'genHTML/chartBoxPlot.js',
							  'genHTML/chartZoom.js',
							  'genHTML/gridism.css',
							  'genHTML/jquery.js',
							  'genHTML/msa.js',
							  'genHTML/scalehd.css',
							  'genHTML/scalehd.js',
							  'genHTML/tablefilter.js',
							  'genHTML/templates/alleletable.html',
							  'genHTML/templates/base.html',
							  'genHTML/templates/fastqc.html',
							  'genHTML/templates/list.html',
							  'genHTML/templates/seqALN.html',
							  'genHTML/templates/seqGTYPE.html',
							  'genHTML/templates/seqqc.html',
							  'genHTML/templates/sequencedata.html',
							  'genHTML/templates/trim.html']},
	include_package_data=True,

	# Executable scripts require an entry point to allow cython to generate
	# executables for the respective target platform. This entry point is akin
	# to launching the script in bash: if __name__ == '__main__' etc..
    entry_points={
        'console_scripts': ['ScaleHD=ScaleHD.sherpa:main',],
    },
)
