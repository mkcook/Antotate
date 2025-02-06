from setuptools import setup
from setuptools import find_packages

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

with open("VERSION.txt", "r") as f:
    version = f.read().rstrip()

with open("requirements.txt", "r") as f:
    requirements = f.readlines()
setup(
    name="antotate",
    version=version,
    author="Margaret Cook",
    author_email="mkcook@uw.edu",
    description="Antotate is a tool allowing for the automatic annotation of Antimony files.",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/mkcook/antotate",
    project_urls={
        "Bug Tracker": "https://github.com/mkcook/antotate/issues",
    },
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    install_requires=requirements,
    scripts=["examples/example.py"],
    package_dir={"": "src"},
    packages=find_packages(where="src"),
    python_requires=">=3.9"
)