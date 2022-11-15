cd NewDepositionData
python3 summarize.py Plans20221106TheHugescan_20221114_v2
cd ..
git status
git add *
git commit -m "Update"
git pull
git push
