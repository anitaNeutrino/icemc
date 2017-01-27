# These aren't the inputs you're looking for...

As I migrated the old `inputs.txt.something` to the *oh-so-fancy* `inputs.something.conf` I found several files that were inconsistent with the "line ordered" variables the old Settings.cc expected.

So at that point I decided to give up recreating all the old files having migrated the major ANITA versions  to the new format.

All the remaining files are left here for archival purposes.

If you really need to use the old file I suggest you copy inputs.anita3.archive.conf as a start and diff the anita3 settings file with the one you want and manually change the variables.
```bash
cp inputs.anita3.archive.conf inputs.something.conf
diff inputs.txt.anita3 inputs.txt.something
```
Since the comments have hardly been altered, the output of the diff should tell you what variables you need to modify.

## Tedious, I know. Good luck :+1:
