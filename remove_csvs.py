#############################################################################################
# Removes all *csv files and allfiles.multi from ddd root directory
#############################################################################################

import subprocess
subprocess.call('rm *csv', shell=True)
subprocess.call('rm allfiles.multi', shell=True)