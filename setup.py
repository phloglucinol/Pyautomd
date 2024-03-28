from os import path
import datetime
from setuptools import setup, find_packages

NAME = 'pyautomd'

readme_file = path.join(path.dirname(path.abspath(__file__)), 'README.md')
try:
    from m2r import parse_from_file

    readme = parse_from_file(readme_file)
except (ImportError, AttributeError) as e:
    with open(readme_file) as f:
        readme = f.read()

today = datetime.date.today().strftime("%b-%d-%Y")
with open(path.join('pyautomd', '_date.py'), 'w') as fp:
    fp.write('date = \'%s\'' % today)


install_requires = [

]

def setup_func(scm=None):
    packages = find_packages(exclude=["example"])

    setup(
        name=NAME,
        use_scm_version=scm,
        setup_requires=['setuptools_scm'],
        author="Runduo Liu",
        author_email="15521120143@163.com",
        description="Python package for automating MD simulations preparation",
        long_description=readme,
        long_description_content_type="text/markdown",
        # url="https://github.com/",
        python_requires="~=3.9",
        packages=packages,
        data_files=[],
        package_data={},
        classifiers=[
            "Programming Language :: Python :: 3.8",
            "License :: OSI Approved :: GNU Lesser General Public License v3 (LGPLv3)",
        ],
        keywords='Pyautomd',
        install_requires=install_requires,
        entry_points={
            'console_scripts': [
            'pyautomd = pyautomd.pyautomd:main'
            ]
        }
    )

try:
    setup_func(scm={'write_to': 'pyautomd/_version.py'})
except:
    setup_func(scm=None)
