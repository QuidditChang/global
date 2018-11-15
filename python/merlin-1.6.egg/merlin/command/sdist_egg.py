"""merlin.command.sdist_egg

Build Python-version-independent source .egg distributions"""


from merlin.command.bdist_egg import bdist_egg
from merlin import get_build_platform, to_filename


class sdist_egg(bdist_egg):


    def finalize_options(self):
        import os
        
        ei_cmd = self.ei_cmd = self.get_finalized_command("egg_info")
        self.egg_info = ei_cmd.egg_info

        if self.bdist_dir is None:
            bdist_base = self.get_finalized_command('bdist').bdist_base
            self.bdist_dir = os.path.join(bdist_base, 'egg')

        if self.plat_name is None:
            self.plat_name = get_build_platform()

        self.set_undefined_options('bdist',('dist_dir', 'dist_dir'))

        if self.egg_output is None:

            # Compute filename of the output egg
            filename = "%s-%s" % (
                to_filename(ei_cmd.egg_name), to_filename(ei_cmd.egg_version)
                )
            platform = self.distribution.has_ext_modules() and self.plat_name
            if platform:
                filename += '-'+self.platform

            self.egg_output = os.path.join(self.dist_dir, filename+'.egg')

        return


    def call_command(self,cmdname,**kw):
        kw.setdefault('compile', False)
        kw.setdefault('optimize', False)
        return bdist_egg.call_command(self, cmdname, **kw)



# end of file
