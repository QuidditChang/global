from merlin.command.easy_install import easy_install
from distutils.util import convert_path
from merlin import Distribution, PathMetadata, normalize_path
from distutils import log
from distutils.errors import *
import sys, os, merlin

class install_deps(easy_install):
    """Install package dependencies."""

    description = "install package dependencies"

    user_options = easy_install.user_options + [
        ("uninstall", "u", "Uninstall this source package"),
    ]

    boolean_options = easy_install.boolean_options + ['uninstall']

    command_consumes_arguments = False  # override base

    def run(self):
        if self.uninstall:
            pass
        else:
            self.install_dependencies()
        self.warn_deprecated_options()

    def initialize_options(self):
        self.uninstall = None
        easy_install.initialize_options(self)










    def finalize_options(self):
        ei = self.get_finalized_command("egg_info")
        if ei.broken_egg_info:
            raise DistutilsError(
            "Please rename %r to %r before using 'install_deps'"
            % (ei.egg_info, ei.broken_egg_info)
            )
        self.args = [ei.egg_name]       
        easy_install.finalize_options(self)
        self.egg_base = ei.egg_base
        self.egg_path = os.path.abspath(ei.egg_base)
        # Make a distribution for the package's source
        self.dist = Distribution(
            normalize_path(self.egg_path),
            PathMetadata(self.egg_path, os.path.abspath(ei.egg_info)),
            project_name = ei.egg_name
        )

    def install_dependencies(self):

        ## NYI: Not sure how to avoid 'egg_info' below... unless the
        ## metadata is written to disk, the dependencies aren't
        ## processed.
        
        # Ensure metadata is up-to-date
        self.run_command('egg_info')
        
        if merlin.bootstrap_install_from:
            self.easy_install(merlin.bootstrap_install_from)
            merlin.bootstrap_install_from = None
        
        # postprocess the installed distro, fixing up .pth, installing scripts,
        # and handling requirements
        self.process_distribution(None, self.dist, not self.no_deps)

    def update_pth(self,dist):
        if dist is not self.dist:
            # Installing a dependency, so fall back to normal behavior
            return easy_install.update_pth(self,dist)

        # no-op
        return

    def install_egg_scripts(self, dist):
        if dist is not self.dist:
            # Installing a dependency, so fall back to normal behavior
            return easy_install.install_egg_scripts(self,dist)

        # no-op
        return
