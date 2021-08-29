from skbuild import setup
from setuptools_scm import get_version

project_name = "GMatElastoPlastic"

setup(
    name = project_name,
    description = "Elasto-plastic material model.",
    long_description = "Elasto-plastic material model.",
    version = get_version(),
    license = "MIT",
    author = "Tom de Geus",
    author_email = "tom@geus.me",
    url = f"https://github.com/tdegeus/{project_name}",
    packages = [f"{project_name}"],
    package_dir = {"": "python"},
    cmake_install_dir = f"python/{project_name}",
    cmake_minimum_required_version = "3.13...3.21",
)
