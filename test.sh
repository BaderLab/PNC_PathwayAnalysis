numperm=$1
covars=$2
for (( i=0; i<=$numperm; i++ ))
do
   echo "$i"
done

if [ $covars = "0" ]
then
   echo "test"
fi
