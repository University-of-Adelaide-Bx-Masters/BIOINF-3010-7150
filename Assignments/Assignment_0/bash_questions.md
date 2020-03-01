# Bash practical assessment task

1. From which directory will `cd ..` not move to the parent directory?

2. Are the paths `~/Practical_1` and `~/Practical_2` relative or absolute paths?

3. Give two ways we could inspect the contents of the `/` directory from your own home directory.

3. The letter `l` and the number `1` are often confused in text, but have different meanings. What is the difference in behaviour of `ls` when run with the `-1` (digit) and `-l` (letter) options? How does `ls -1` (digit) differ from `ls` without options?

4. If we wanted to hide the group names in the long listing format, which extra options would we need set when searching our home directory?

5. Try accessing the documentation for the command `man` all the ways you can think of. Was there a difference in the output depending on how we asked to view the documentation? Could you access the documentation for the `ls` command all three ways?

6. Complete the table.

	| **Command** | **Description of function**   |
	|:----------- |:----------------------------- |
	| `mv`        |                               |
	| `cp`        |                               |
	| `rm`        |                               |
	| `mkdir`     |                               |
	| `cat`       |                               |
	| `less`      |                               |
	| `wc`        |                               |
	| `head`      |                               |
	| `tail`      |                               |
	| `echo`      |                               |
	| `cut`       |                               |
	| `sort`      |                               |
	| `uniq`      |                               |
	| `wget`      |                               |
	| `gunzip`    |                               |

7. What is "stdin" an abbreviation of?

8. Why does `echo ~` output `/home/student` when you execute it? What would it output if a user with a different home directory executed it? What happens when you execute `echo ~/*`?

9. How many features are contained in the GCF\_000182855.2\_ASM18285v1\_genomic.gff file?
Why is `wc -l GCF_000182855.2_ASM18285v1_genomic.gff` not a correct way to determine the number of features?

10. What does the `'^>'` mean in the grep command where Drosophila\_melanogaster.BDGP6.ncrna.fa is being searched?

11. Why doesn't `grep -c 'gene' GCF_000182855.2_ASM18285v1_genomic.gff` print the number of genes described in the GCF\_000182855.2\_ASM18285v1\_genomic.gff file?

12. In the first script that you used in the practical, there are two variables. What are their names?

13. Before you executed the `chmod` command on the wellDone.sh file, what were the permissions? What did these permissions allow you to do as a user, as members of your user group, as another user not in your user group?

14. What did the 4 digit do to the permissions of wellDone.sh with the command `chmod 774 wellDone.sh`?

15. Add [meaningful comments](http://doc.cat-v.org/bell_labs/pikestyle) to the following small script (I would suggest using more than Rob does for this exercise - use your judgement):

	```
	#!/bin/bash

	FILES=$(ls)

	COUNT=0
	for f in ${FILES}; do
		((COUNT++))
		ln=$(wc -l ${f} | cut -f 1 -d ' ')
		echo "File number ${COUNT} (${f}) has ${ln} lines"
	done
	```

	Note that you may not have comments above the shebang line; `#!` must be the first two characters in the file for it to function correctly.
