

from merlin import Command
#from merlin.command.easy_install import easy_install


class egg_flags(Command):

    """Command to get egg-related build flags."""

    description = "pring egg-related build flags"

    user_options = [
    ]


    def initialize_options(self): pass


    def finalize_options(self):
        import os
        from merlin import Distribution, PathMetadata, normalize_path

        ei = self.get_finalized_command("egg_info")
        if ei.broken_egg_info:
            raise DistutilsError(
            "Please rename %r to %r before using 'egg_flags'"
            % (ei.egg_info, ei.broken_egg_info)
            )
        self.args = [ei.egg_name]       
        #easy_install.finalize_options(self)
        self.egg_base = ei.egg_base
        self.egg_path = os.path.abspath(ei.egg_base)
        # Make a distribution for the package's source
        self.dist = Distribution(
            normalize_path(self.egg_path),
            PathMetadata(self.egg_path, os.path.abspath(ei.egg_info)),
            project_name = ei.egg_name
        )


    def run(self):

        from merlin import WorkingSet

        import sys
        from ConfigParser import ConfigParser, NoOptionError
        from StringIO import StringIO

        flags = dict(
            CFLAGS = [],
            CPPFLAGS = [],
            LDFLAGS = [],
            LIBS = [],
            PYXFLAGS = [],
            )

        requirements = [self.dist.as_requirement()]
        working_set = WorkingSet()

        deps = working_set.resolve(requirements)
        deps.reverse()
        dependencies = []
        processed = {}
        for dist in deps:
            if dist in processed:
                continue
            dependencies.insert(0, dist)
            processed[dist] = True
        for dist in dependencies:
            if dist.has_metadata('config.cfg'):
                parser = ConfigParser({'location': dist.location})
                config = dist.get_metadata('config.cfg')
                fp = StringIO(config)
                parser.readfp(fp, 'config.cfg')
                for k,v in flags.iteritems():
                    try:
                        v.append(parser.get('flags', k))
                    except NoOptionError:
                        pass

        stream = open("egg-flags.sh", "w")
        for k,v in flags.iteritems():
            print >> stream, 'PYTHON_EGG_%s="%s"' % (k, ' '.join(v))
        stream.close()

        return


# end of file
