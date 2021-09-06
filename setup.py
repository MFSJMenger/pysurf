#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""The setup script."""

from setuptools import setup, find_packages

with open('README.rst') as readme_file:
    readme = readme_file.read()

with open('HISTORY.rst') as history_file:
    history = history_file.read()

requirements = [ ]

setup_requirements = ['pytest-runner', ]

test_requirements = ['pytest>=3', ]

setup(
    author="Maximilian Menger, Johannes Ehrmaier",
    author_email='menger.maximilian@gmail.com',
    python_requires='>=2.7, !=3.0.*, !=3.1.*, !=3.2.*, !=3.3.*, !=3.4.*, !=3.5.*',
    classifiers=[
        'Development Status :: 2 - Pre-Alpha',
        'Intended Audience :: Developers',
        'License :: OSI Approved :: Apache License 2.0',
        'Natural Language :: English',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
    ],
    description="Python Surface Hopping Code",
    entry_points={
        'console_scripts': [],
    },
    install_requires=requirements,
    license="Apache License 2.0",
    long_description=readme + '\n\n' + history,
    include_package_data=True,
    keywords='pysurf',
    name='pysurf',
    packages=find_packages(include=['pysurf', 'pysurf.*']),
    setup_requires=setup_requirements,
    test_suite='tests',
    tests_require=test_requirements,
    url='https://github.com/mfsjmenger/pysurf',
    version='0.1.0',
    zip_safe=False,
)
