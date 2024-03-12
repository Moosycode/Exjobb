time_seconds = [
    10026, 10034, 10041, 10049, 10057, 10065, 10072, 10080, 10088, 10095,
    10103, 10111, 10119, 10126, 10134, 10142, 10150, 10157, 10165, 10173,
    10181, 10188, 10196, 10204, 10212, 10219, 10227, 10235, 10242, 10250,
    10258, 10266, 10273, 10281, 10289, 10297, 10304, 10312, 10320, 10328,
    10335, 10343, 10351, 10359, 10366, 10374, 10382, 10389, 10397, 10405,
    10413, 10420, 10428, 10436, 10444, 10451, 10459, 10467, 10475, 10482,
    10490, 10498, 10506, 10513, 10521, 10529, 10536, 10544, 10552, 10560,
    10567, 10575, 10583, 10591, 10598, 10606, 10614, 10622, 10629
]
current_minute = -1  # Initialize current minute
for seconds in time_seconds:
    minute = seconds // 60  # Calculate current minute
    if minute != current_minute:  # Check if minute has changed
        print("A minute has passed at", seconds, "seconds.")
        current_minute = minute  # Update current minute