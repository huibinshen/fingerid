

from setuptools import setup, find_packages

config = {
    'description':'fingerid-package',
    'author':'Huibin Shen',
    'url':'project https://sourceforge.net/projects/fingerid/',
    'author_email':'huibin.shen@aalto.fi',
    'version':'2.0',
    'install_requires':['nose'],
    'packages':find_packages(),
    'name':'fingerid',
}

setup(**config)
