import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="SPOCK", # Replace with your own username
    version="0.0.1",
    author="Elsa Ducrot",
    author_email="educrot@uliege.be",
    long_description = read('README.rst'),
    install_requires = ['pandas','numpy==1.19.0','astroplan','astropy','matplotlib','datetime','pyaml','docx','plotly',
                    'gspread, oauth2client, astroplan, alive_progress, paramiko , chart_studio, python-docs, markdown'],
    url = 'https://github.com/educrot/SPOCK/',
    packages=setuptools.find_packages(),
    python_requires='>=3.6',
)
