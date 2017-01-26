# These aren't the inputs you're looking for...

As I migrated the old inputs.txt.something to the oh so fancy inputs.something.conf I found several files that were inconsistent with the line ordered variables the old Settings.cc expected.
So I decided to give up at that point having migrated a few of the major ANITA version inputs to the new format.
All the remaining files are left here for archival purposes.
If you really need to use the old file I suggest:
	- you copy inputs.anita3.archive.conf
	- diff the file you want to migrate to the new format with inputs.txt.anita3
	- Since the comments have hardly been altered, the output of the diff should tell you what variables you need to modify.
Tedious, I know. Good luck.
