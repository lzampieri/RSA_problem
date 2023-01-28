cd NewDepositionData
python3 summarize.py Plans20230110OnlySquaredonlyQ02scan_20230109_v2
cd ..
git status
git add *
git commit -m "Update"
git pull
git push
