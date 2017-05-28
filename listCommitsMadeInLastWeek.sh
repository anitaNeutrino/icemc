# Script to print the commits to the software framework made in the last week.

theOriginalPwd=$(pwd);

for component in ./ ../anitaBuildTool/RandomMacros/createIceMCinputs/; do
    cd $component;
    echo $component | sed s,"components/",,;

    # First get a list of all hashes from all branches from the last 7 days    
    theHashes=$(git log --branches --all --remotes=* --since=7.days.ago --pretty=format:"%h")

    for theHash in ${theHashes}; do
	# List all branches with a tip that is a descendent of that commit (include merges I guess)
	# -r since we only care about remotely tracked branches (and I branch for almost every new feature locally)
	# The two sed commands remove the HEAD branch pointer and the origin prefix from the remote branch
	theBranches=$(git branch -r --contains ${theHash} | sed s,"origin/HEAD -> origin/master",, | sed s,"origin/",,)

	# Print all the branches
	echo -n ${theBranches} ""
	git --no-pager show ${theHash} --pretty=format:"%h%x09%an%x09%s" --summary
    done
    echo "";
    echo "";
    cd $theOriginalPwd
done;    
