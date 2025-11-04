echo $1
convert $1.png -channel RGB -negate $1_light.png
convert $1_light.png -fuzz 5% -fill '#131416' -opaque black $1_dark.png
rm $1_light.png
