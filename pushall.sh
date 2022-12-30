cd NewDepositionData
python3 summarize.py Plans20221212OnlyDimersonlyQ02scan_20221212_v2
cd ..
git status
git add *
git commit -m "Update"
git pull
git push
