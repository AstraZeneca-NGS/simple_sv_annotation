#!/usr/bin/env python
from setuptools import setup

package_name = 'simple_sv_annotation'

setup(
    name=package_name,
    version='1.0',
    description='Prioritizing SV variants',
    keywords='bioinformatics',
    license='GPLv3',
    include_package_data=True,
    zip_safe=False,
    install_requires=['pyvcf'],
    scripts=['simple_sv_annotation.py'],
    classifiers=[
        'Environment :: Console',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
        'Natural Language :: English',
        'Operating System :: MacOS :: MacOS X',
        'Operating System :: POSIX',
        'Operating System :: Unix',
        'Programming Language :: Python',
        'Programming Language :: JavaScript',
        'Topic :: Scientific/Engineering',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
    ],
)