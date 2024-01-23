import shutil

class Executable_checker:
    def __init__(self, ):
        pass

    def is_executable_available(self, program_name):
        """Checks if the program is available in the system PATH and executable."""
        return shutil.which(program_name) is not None


