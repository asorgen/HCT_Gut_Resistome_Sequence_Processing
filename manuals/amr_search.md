# AMR detection

**Description:**

## Module overview


## Setup

See [setup](manuals/setup.md) for requirements. 


## Run module

Assuming you have this github cloned, you will need to enter the root cloned directory and enter the following command:

**Short-read analysis**
```shell
nohup sh ./SHORT/short_pipeline.sh amr_search > SHORT/short_pipeline.out 2>&1 &
```

**Long-read analysis**
```shell
nohup sh ./LONG/long_pipeline.sh amr_search > LONG/long_pipeline.out 2>&1 &
```

**Hybrid assembly analysis**
```shell
nohup sh ./HYBRID/hybrid_pipeline.sh amr_search > HYBRID/hybrid_pipeline.out 2>&1 &
```


## Output

**This module creates the following output in your SHORT, LONG, OR HYBRID folders:**

