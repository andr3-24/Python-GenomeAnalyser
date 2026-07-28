from pathlib import Path
from datetime import datetime

#get bin
BASE_DIR = Path(__file__).resolve().parents[1] 


#Get log folder
log_file = BASE_DIR / "Bin" / "Log" / "log.txt"


def dateTime():
    return datetime.now().strftime("%d-%m-%Y_%H-%M-%S")
        
        
def initLog():
    # Clear or create log file
    log_file.write_text(
        f"New log started {dateTime()}:\n",
        encoding="utf-8"
    )


def log(*args):
    with open(log_file, "a", encoding="utf-8") as file:
        print(dateTime(), ">>>", *args, file=file)
      
        