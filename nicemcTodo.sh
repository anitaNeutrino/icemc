echo "***************************"
echo "****THINGS MARKED TODO*****"
echo "***************************"
grep --color=always -n "@todo" src/* include/*
echo "***************************"
echo ""
echo ""
echo ""
echo "***************************"
echo "****THINGS MARKED DEBUG****"
echo "***************************"
grep --color=always -n "@attention" src/* include/*
echo "***************************"
