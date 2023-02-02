for f in *.zip; do zip -d "$f" "__MACOSX/*"; done
for f in *.zip; do zip -d "$f" "*/.DS_Store"; done
