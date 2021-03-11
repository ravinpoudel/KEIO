"""setup.py: python package setup for KEIO

"""

from setuptools import setup

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
    keywords='KIEIO,
    url='https://ravinpoudel.github.io',
    test_suite='pytest',
    author='Ravin Poudel',
    author_email='rp3448@ufl.edy',
    install_requires=['biopython>=1.70', 'pybedtools>=0.8.0', 'nmslib>=2.0.4', 'pandas>=1.0.0'],
    python_requires='>=3.6',
    tests_require=['pytest'],
    include_package_data=True,
    entry_points={'console_scripts':['keio = keio.cli:main',]},
    zip_safe=False)
