import sys,commands

class StageError(Exception):
    def __init__(self, value):
        self.value = value
    def __str__(self):
        return repr(self.value)

def checkStage(stage):
    #print("Checking:",stage)
    cmd = commands.getoutput('./modmul stage'+str(stage)+' < stage'+str(stage)+'.input')
    cmd = cmd.split()

    if (len(cmd)==0):
        raise StageError('This stage did not return any output')

    for idx, item in enumerate(cmd):
        cmd[idx] = cmd[idx]
        cmd[idx] = cmd[idx].upper()

    stagename = 'stage'+str(stage)+'.output'
    correct = open(stagename, 'r')
    correct = correct.read()
    correct = correct.split()

    for idx, item in enumerate(correct):
        if (correct[idx] != cmd[idx]):
            print('Error in stage %d, case %d'%(stage,idx))
            print('Test: ',cmd[idx])
            print('Correct: ',correct[idx])
            return False 
    return True

for i in range(1,5):
    try:
        if checkStage(i):
            print('Stage '+str(i)+': PASS')
    except StageError as e:
        print('Stage '+str(i)+': '+e.value)