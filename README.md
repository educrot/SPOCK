
![Test Image 1](logo_SPOCK_2.png)


# SPOCK

**`SPOCK`** is a Python library for dealing with the planification of SPECULOOS targets observations

*Schedule targets on several criteria:*
*  Visibility of the target
*  Priority (metric from JWST)
*  number of hours already performed

## Installation

Use the package manager [pip]() to install SPOCK.

```bash
pip install SPOCK
```

## Usage

Create your `'input_file.csv'` file in the following format:

--- 
    date_range: 
      - "2019-09-23 15:00:00"
      - "2019-09-26 15:00:00"
    duration_segments: 20
    nb_segments: 3
    observatories: 
      - Saint-Ex
      - SNO
    strategy: "continuous"
    target_list: target_prio_50.txt
---

Then, open a python script or the [SPOCK jupyter notebook]() and run:

```python
import SPOCK.long_term_scheduler as SPOCKLS
import SPOCK.plots_scheduler as SPOCKplots

schedule = SPOCKLS.schedule()
schedule.load_parameters('input_file.csv')
schedule.make_schedule()
schedule.make_plan_file('input_file.csv')
```

To plot the schedule you  have generated, execute the following command:

```
fig = SPOCKplots.gantt_chart('plan.txt',schedule.observatory)
```

Example of output image you will obtain:


![Test Image 1](schedule_example.png)


## Contributing
Pull requests are welcome. For major changes, please open an issue first to discuss what you would like to change.

Please make sure to update tests as appropriate.

## License

<span style=“color:red;”> text </span>