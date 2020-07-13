# -*- coding: utf-8 -*-

from setuptools import setup, find_packages

with open('README.md') as f:
    readme = f.read()

with open('LICENSE.md') as f:
    license = f.read()

setup(
    name='psaap3reac0d',
    version='0.0.1',
    description='Reduced Order Modeling for the PSAAP3 project.',
    long_description=readme,
    author='Quentin Douasbin',
    author_email='quentin.douasbin@gmail.com',
    url='',
    license=license,
    packages=find_packages(exclude=('tests', 'docs'))
)
