"""setup.py: python package setup for KEIO

"""

from setuptools import setup

requirements = [
    # package requirements go here
    'biopython>=1.79',
    'numpy >=1.11',
    'pybedtools>=0.8.2',
    'nmslib>=2.0.6',
    'pandas>=1.0.0',
    'pytest>=4.6',
    'textdistance>=4.2.1',
    'pytest-cov'
]

setup(
    name='KEIO',
    packages=['KEIO'],
    license='CC0 1.0 Universal (CC0 1.0) Public Domain Dedication',
    description='KEIO: A python software to process illumina reads for keio-collection type project',
    long_description=open('README.md').read(),
    classifiers=['Topic :: Scientific/Engineering :: Bio-Informatics',
                 'Programming Language :: Python :: 3.7',
                 'Programming Language :: Python :: 3.8',
                 'Development Status :: 3 - Alpha'],
    keywords='KEIO',
    url='https://ravinpoudel.github.io',
    test_suite='pytest',
    author='Ravin Poudel',
    author_email='rp3448@ufl.edu',
    install_requires=requirements,
    python_requires='>=3.6',
    tests_require=['pytest'],
    include_package_data=True,
    entry_points={'console_scripts':['keio = keio.cli:main',]},
    zip_safe=False)
