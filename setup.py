import os
from setuptools import setup, find_packages

# Utility function to read the README file.
def read(fname):
    return open(os.path.join(os.path.dirname(__file__), fname)).read()

setup(
    name = 'SPOCK',
    version = '0.0.0',
    author = 'Elsa Ducrot',
    author_email = 'educrot@uliege.be',
    description = ('Speculoos Observatory SChedule maKer'),
    keywords = '',
    url = 'https://github.com/educrot/SPOCK/',
    packages = find_packages(),
    long_description = read('README.rst'),
    install_requires = ['pandas','numpy==1.19.0','astroplan','astropy','matplotlib','datetime','pyaml','docx','plotly'],
    python_requires='>=3.6',
)