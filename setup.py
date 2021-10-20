import os
from setuptools import setup, find_packages

# Utility function to read the README file.
def read(fname):
    return open(os.path.join(os.path.dirname(__file__), fname)).read()

setup(
    name = 'SPOCK',
    version = '0.0.1',
    author = 'Elsa Ducrot',
    author_email = 'educrot@uliege.be',
    description = ('Speculoos Observatory SChedule maKer for chilean night on SPECULOOS South Observatory'),
    keywords = '',
    url = 'https://github.com/educrot/SPOCK/',
    packages = find_packages(),
    long_description = read('README.rst'),
    python_requires='>=3.6',
    install_requires =['pandas','tqdm','numpy','astroplan','astropy','matplotlib','datetime','pyaml',
                       'plotly','colorama','gspread', 'oauth2client', 'astroplan', 'alive_progress', 'paramiko',
                       'requests','chart_studio', 'markdown','python-docx','bs4','ipywidgets==7.6.3'],
)