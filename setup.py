# Always prefer setuptools over distutils
from setuptools import setup, find_packages
from codecs import open
from os import path

here = path.abspath(path.dirname(__file__))

# Get the long description from the README file
with open(path.join(here, 'README.md'), encoding='utf-8') as f:
    long_description = f.read()

setup(
    name='basematic',
    version='0.1.1',
    description='Tools For Bioinformatics',
    long_description=long_description,
    long_description_content_type='text/markdown',
    url='https://gene.pku.edu.cn',
    author='Xiannian Zhang',
    author_email='friedpine@gmail.com',

    classifiers = [
        'Development Status :: 3 - Alpha',
        'Intended Audience :: Developers',
        'Topic :: Software Development :: Build Tools',
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python :: 3.6',
    ],

    keywords = 'sample setuptools development',

    packages = find_packages(include=['basematic*']),

    install_requires=[
        'pandas',
        'numpy',
        'matplotlib',
        'jinja2',
        'click',
        'configparser'
    ],

    extras_require={
        'dev': ['check-manifest'],
        'test': ['coverage'],
    },

    entry_points={
        'console_scripts': [
            'basematic=basematic.cmd:cli',
            'basematic-manage=basematic.mgt:cli',
            'basematic-fastq=basematic.bio.fastq:cli',
            'basematic-CNV=basematic.bio.cnv:cli',
            'basematic-BAM=basematic.bio.bam:cli',
            'basematic-SNV=basematic.bio.snv._cmd:cli',
            'basematic-RNA=basematic.bio.rnaseq:cli',
        ]
    }
)