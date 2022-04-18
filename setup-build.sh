#! /usr/bin/env bash

data="$(cat)"
alphabet="$(echo "$data" | head -n 1)"
prec="$(echo "$data" | head -n 2 | tail -n 1)"
words="$(echo "$data" | tail -n +3)"

n=0
for i in $alphabet; do
        n=$((n+1))
        eval "a$i=$n"
done
maxl=$(
echo "$words" | while read ww; do
        echo "$ww" | wc -w
done | sort -n | tail -n 1
)
echo "$alphabet" | wc -w
echo "$alphabet"
echo "$((maxl-1))"
echo "$words" | wc -l
echo "$words" | while read ww; do
        echo -n "$(echo "$ww" | wc -w) "
        for i in $ww; do
                eval "echo -n \$a$i ''"
        done
        echo
done
echo $prec
