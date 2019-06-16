import os
import re
import sys
import platform
import subprocess

from setuptools import setup, Extension, find_packages
from setuptools.command.build_ext import build_ext
from distutils.version import LooseVersion

pipe = subprocess.Popen('git describe --tags --always',
                        stdout=subprocess.PIPE, shell=True)
git = pipe.stdout.read().decode("utf-8").rstrip()
git_release = git.lstrip('v')
git_version = '.'.join(git_release.split('-')[0:2])

cpus = os.cpu_count()

here = os.path.abspath(os.path.dirname(__file__))

with open(os.path.join(here, 'README.md'), encoding='utf-8') as f:
    long_description = f.read()

with open(os.path.join(here, 'requirements.txt')) as f:
    requirements = f.read().split()


class CMakeExtension(Extension):
    def __init__(self, name, sourcedir=''):
        Extension.__init__(self, name, sources=[])
        self.sourcedir = os.path.abspath(sourcedir)


class CMakeBuild(build_ext):
    def run(self):
        try:
            out = subprocess.check_output(['cmake', '--version'])
        except OSError:
            raise RuntimeError("CMake must be installed to build the following"
                               "extensions: " +
                               ", ".join(e.name for e in self.extensions))

        if platform.system() == "Windows":
            cmake_version = LooseVersion(
                re.search(r'version\s*([\d.]+)', out.decode()).group(1))
            if cmake_version < '3.1.0':
                raise RuntimeError("CMake >= 3.1.0 is required on Windows")

        for ext in self.extensions:
            self.build_extension(ext)

    def build_extension(self, ext):
        extdir = os.path.abspath(
            os.path.dirname(self.get_ext_fullpath(ext.name)))
        cmake_args = ['-DCMAKE_LIBRARY_OUTPUT_DIRECTORY=' + extdir,
                      '-DPYTHON_EXECUTABLE=' + sys.executable]

        cfg = 'Debug' if self.debug else 'Release'
        build_args = ['--config', cfg]

        if platform.system() == "Windows":
            cmake_args += ['-DCMAKE_LIBRARY_OUTPUT_DIRECTORY_{}={}'
                           .format(cfg.upper(), extdir)]
            if sys.maxsize > 2**32:
                cmake_args += ['-A', 'x64']
            build_args += ['--', '/m']
        else:
            cmake_args += ['-DCMAKE_BUILD_TYPE=' + cfg]
            cmake_args += ['-DDEBUG_OUTPUT=ON']
            # cmake_args += ['-DCMAKE_VERBOSE_MAKEFILE:BOOL=ON']
            build_args += ['--target', '_pycentroids']
            build_args += ['--', '-j{}'.format(cpus)]

        env = os.environ.copy()
        env['CXXFLAGS'] = '{} -DPY_VERSION_INFO=\\"{}\\"' .format(
            env.get('CXXFLAGS', ''),
            self.distribution.get_version())
        if not os.path.exists(self.build_temp):
            os.makedirs(self.build_temp)
        subprocess.check_call(['cmake', ext.sourcedir] + cmake_args,
                              cwd=self.build_temp, env=env)
        subprocess.check_call(['cmake', '--build', '.'] + build_args,
                              cwd=self.build_temp)


setup(
    name='pycentroids',
    author='Stuart B. Wilkins',
    author_email='swilkins@bnl.gov',
    description='Centroiding algorithms for CCD Single Photon Counting',
    long_description=long_description,
    long_description_content_type='text/markdown',
    license='BSD (3-clause)',
    url='https://github.com/NSLS-II/centroids',
    packages=find_packages(),
    setup_requires=["flake8"],
    install_requires=requirements,
    ext_modules=[CMakeExtension('_pycentroids')],
    cmdclass=dict(build_ext=CMakeBuild),
    zip_safe=False,
    version=git_version
)
