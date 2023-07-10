import setuptools
setuptools.setup(name='subproptools',
version='0.0.2',
description='Extract group properties from .sum file and rotate substituent geometries',
url='https://github.com/kmlefran/subproptools',
author='Kevin M. Lefrancois-Gagnon',
install_requires=['numpy','pandas'],
author_email='kgagnon@lakeheadu.ca',
packages=setuptools.find_packages(),
zip_safe=False)