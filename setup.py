

from setuptools import setup, find_packages

config = {
    'description':'fingerid-package',
    'author':'Huibin Shen',
    'url':'project https://github.com/icdishb/fingerid',
    'author_email':'huibin.shen@aalto.fi',
    'version':'1.4',
    'install_requires':['nose'],
    'packages':find_packages(),
    'name':'fingerid',
}

setup(**config)
