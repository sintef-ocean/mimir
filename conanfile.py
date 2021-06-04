from conans import ConanFile, CMake, tools
import os
# import json
# import pprint


class MimirConan(ConanFile):
    name = "mimir"
    license = "Apache-2.0"
    url = "https://github.com/sintef-ocean/mimir"
    author = "SINTEF Ocean"
    homepage = "https://sintef-ocean.github.io/mimir"
    description = \
        "Mimir is an path planning application for purse seining"
    topics = ("QML", "DDS", "OpenSplice", "Purse seining", "Nonlinear programming")
    exports = "version.txt, LICENSE"
    settings = "os", "compiler", "build_type", "arch"
    keep_imports = True
    generators = (
      "virtualrunenv",
      "virtualenv",
      "cmake",
      "cmake_paths",
      "cmake_find_package",
      "json"
    )
    options = {
      "fPIC": [True, False],
      "with_doc": [True, False],
      "with_gnuplot": [True, False],
      "with_CICD": [True, False]
    }
    default_options = (
      "fPIC=True",
      "with_doc=False",
      "with_gnuplot=False",
      "with_CICD=False"
    )

    build_subfolder = "build_subfolder"
    _cmake = None

    def requirements(self):
        self.requires("casadi/[>=3.5.5]@sintef/stable")
        self.requires("opensplice-ce/[>=6.9]@sintef/stable")
        self.requires("RatatoskIDL/0.2.0@sintef/stable")
        self.requires("yaml-cpp/[>=0.6.0]")
        self.requires("boost/1.69.0")

        if self.options.with_gnuplot:
            self.requires("gnuplot-iostream/2020.06.20@sintef/stable")

    def build_requirements(self):
        # Internal Continuous deployment helper scripts
        if self.options.with_CICD:
            self.build_requires("kluster-scripts/[>=0.2.0]@kluster/stable")

        if tools.os_info.is_linux:
            installer = tools.SystemPackageTool()
            installer.install("patchelf")

    def configure(self):

        self.options["boost"].shared = True
        self.options["yaml-cpp"].shared = True
        self.options["casadi"].thread = True
        self.options["openblas"].dynamic_arch = True
        # self.options["openblas"].use_thread = True
        self.options["casadi"].with_common = True

        if self.options.with_gnuplot:
            self.options["gnuplot-iostream"].with_boost = False
            self.options["zstd"].shared = True

    def export_sources(self):

        self.copy("*", src="cmake", dst="cmake")
        self.copy("*", src="data", dst="data")
        self.copy("*", src="docs", dst="docs")
        self.copy("*", src="examples", dst="examples")
        self.copy("*", src="external", dst="external")
        self.copy("*", src="include", dst="include")
        self.copy("*", src="src", dst="src")
        self.copy("*", src="test_package", dst="test_package")
        self.copy("*", src="tools", dst="tools")
        self.copy("CMakeLists.txt")
        self.copy("version.txt")
        self.copy("LICENSE")
        self.copy("README.org")

    def set_version(self):
        self.version = tools.load(os.path.join(self.recipe_folder,
                                               "version.txt")).strip()

    def _configure_cmake(self):
        if self._cmake:
            return self._cmake

        self._cmake = CMake(self)
        self._cmake.definitions["MIMIR_WITH_GNUPLOT"] = self.options.with_gnuplot
        self._cmake.definitions["WITH_DOC"] = self.options.with_doc

        if self.settings.os != "Windows":
            self._cmake.definitions["CMAKE_POSITION_INDEPENDENT_CODE"] = self.options.fPIC

        self._cmake.configure()
        return self._cmake

    def build(self):
        cmake = self._configure_cmake()
        cmake.build()

        if self.options.with_doc:
            cmake.build(target='doc')

    def package(self):
        cmake = self._configure_cmake()
        cmake.build(target='package_it')

        self.copy("{}*.tar.gz".format(self.name), dst=".", keep_path=False)
        self.copy("{}*.exe".format(self.name), dst=".", keep_path=False)
        self.copy("{}*.deb".format(self.name), dst=".", keep_path=False)

    def imports(self):

        self.copy("*", src="licenses",
                  dst="bundle/share/mimir/licenses",
                  folder=True)

        copy_patterns = ["*.so*", "*.dll"]
        copy_sources = ["@bindirs", "@libdirs"]
        boost_libs = [
            "*program_options*",
            "*log*",
            "*filesystem*",
            "*thread*",
            "*date_time*"
        ]
        if self.settings.os == "Windows":
            bundle_dest = "bundle/bin"
        else:
            bundle_dest = "bundle/lib"
        if self.options.with_gnuplot:
            boost_libs.append("*iostreams*")
            boost_libs.append("*system*")

        for copy_pattern in copy_patterns:
            for copy_source in copy_sources:
                self.copy(copy_pattern,
                    src=copy_source,
                    dst=bundle_dest,
                    excludes=["*boost*"],
                    keep_path=False)

        for b_lib in boost_libs:
            for copy_source in copy_sources:
                self.copy(b_lib,
                        src=copy_source,
                        root_package="boost",
                        dst=bundle_dest,
                        excludes=["*.lib"],
                        keep_path=False)

        # Copy default OpenSplice configuration files
        self.copy("*ospl.xml",
                  dst="bundle",
                  root_package="opensplice-ce",
                  keep_path=True)
        self.copy("*ospl_metaconfig.xml",
                  dst="bundle",
                  root_package="opensplice-ce",
                  keep_path=True)

    def package_id(self):
        del self.info.settings.compiler
        # other stuff that are not relevant for runtime package
