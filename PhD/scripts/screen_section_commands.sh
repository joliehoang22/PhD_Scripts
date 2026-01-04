create: screen -S screen_section_name
detach: Ctrl + a (release) d
list all screen section: screen -ls
enter a screen section: screen -r screen_section_name
delete a screen section within a section: Ctrl + a + k
help: Ctrl + a + ?
delete a screen section outside: screen -S 1899356.brainsegfounder -X quit

lshosts #list all the hosts 


screen -S 3874189.3738434.shuffling -X quit
screen -S 3744941.shuffling2 -X quit
screen -S 3738434.shuffling -X quit

#to view a file in the terminal
cat mini_fold.json #cat = “concatenate and print” → shows the file contents directly.
less mini_fold.json # scroll through the file
head mini_fold.json 
tail mini_fold.json 




