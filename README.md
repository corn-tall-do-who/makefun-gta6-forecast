[![make.fun · gta6-twitch-ccu-forecast-01](https://make.fun/badge/gta6-twitch-ccu-forecast-01.svg)](https://make.fun/bounty/gta6-twitch-ccu-forecast-01?utm_source=github&utm_medium=readme&utm_campaign=gta6-twitch-ccu-forecast-01)
[![rank](https://make.fun/badge/gta6-twitch-ccu-forecast-01/corn-tall-do-who/makefun-gta6-forecast/rank.svg)](https://make.fun/bounty/gta6-twitch-ccu-forecast-01/leaderboard?utm_source=github&utm_medium=readme&utm_campaign=gta6-twitch-ccu-forecast-01)
[![reward](https://make.fun/badge/gta6-twitch-ccu-forecast-01/corn-tall-do-who/makefun-gta6-forecast/reward.svg)](https://make.fun/bounty/gta6-twitch-ccu-forecast-01?utm_source=github&utm_medium=readme&utm_campaign=gta6-twitch-ccu-forecast-01)

# GTA VI Twitch CCU forecaster — variant B (corn-tall-do-who)

Forecasts how many people will be watching Grand Theft Auto VI live on
Twitch in each hour of its release day (hours 2–24, UTC), as five quantiles
an hour, from the first hour's viewer count alone.

## The model

The forecast is anchored on the first hour's print: every later hour is a
multiple of it. The multiple's centre — the median — is a piecewise-linear
function of the hour in log space, with three segments (early hours, middle
of the day, late hours), each with its own slope. The slopes are fitted at
run time from the past big launches the case file itself carries: each
comparable's release-day curve is divided by its own hour-1 count and
logged, and only launches whose hour-1 count is within a factor of five of
this case's hour-1 count are used as peers. A launch that had barely
started at 00:00 UTC behaves nothing like one that was already at full
launch level, and the anchor's size is what tells them apart; when too few
peers qualify, the band widens fivefold at a time before falling back to
all usable comparables.

The spread widens with the hour. Around the fitted median, the peer
launches' residuals are pooled over a five-hour window centred on the hour
being forecast, and each level's offset is that window's empirical
quantile of the residuals — so the outer bands are honest about how far
peer launches actually strayed, including the long upper tail of launches
that unlock later in the day, rather than assuming a symmetric bell shape.

## Local validation

Held out, one launch at a time, over the 37 comparable launches in the
practice case (each launch forecast from its own hour-1 print with itself
removed from the peer pool), the forecaster's mean pinball loss is 0.398
against 1.142 for the flat baseline that holds the hour-1 print — about
65% lower — and 0.362 against 0.894 on the 34 launches whose hour-1 count
is at least 1,000. It beats the flat baseline on 36 of the 37 launches.
Two runs on the same input are byte-identical, and all four practice cases
run in about a second in total, far inside the 8-second-per-case cap.

## Running it

```sh
MAKEFUN_INPUT_DIR=input MAKEFUN_OUTPUT_DIR=out python3 main.py
```

Python 3.13 standard library only, no network, deterministic. Reads every
`*.in` in the input directory and writes one `*.out` per case to the output
directory: `{"hours": [{"hour": h, "q": [q05, q25, q50, q75, q95]}]}` for
hours 2–24.

## License

MIT — see `LICENSE`.
