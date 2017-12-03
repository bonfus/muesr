import sys
from muesr.core.nprint import nprint

try:
    import readline
except:
    nprint ("readline not present, using standard python input functions.\n",'warn')


def ninput(message,parser = None):
    """
    Nice input function
    """
    try:
        choice = ""
        if sys.version_info >= (3,0):
            choice = input('\t ' + message )
        else:
            choice = raw_input('\t ' + message )
        #readline.remove_history_item(readline.get_current_history_length()-1)
        if parser:
            choice = parser(choice)
        return choice
    except KeyboardInterrupt:
        nprint ("Use crtl+D to stop input!\n",'warn')
        pass

def ninput_mt(message, parser = None, emessage = 'Parser could not parse your input or invalid input.'):
    """
    Nice input multiple trials function
    """    

    
    while True:
        try:
            ui = ninput(message)
            if parser:
                ui = parser(ui)
            return ui
        except ValueError:  # tries again if parser fails
            nprint(emessage,'warn')
            

