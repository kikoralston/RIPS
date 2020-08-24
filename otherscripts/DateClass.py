# Michael Craig, 6 June 2016

# Date class


# INPUTS: year, month, and day as strings or integers;
# hour as 1-24 int or HH:HH str;
# timezone (available: UTC, EST).
class Date(object):
    daysPerMonth = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]

    def __init__(self, year, month, day, hour, timezone):
        self.year = int(year)
        if self.year < 100: self.year += 2000  # entered as 06
        self.month = int(month)
        self.day = int(day)
        if type(hour) == int:
            self.hour = hour
        else:
            if ':' not in hour:
                self.hour = int(hour)
            else:
                self.hour = int(hour[:hour.index(':')]) + 1  # add 1 to convert to 1-24
        self.timezone = timezone

    def __eq__(self, other):
        if type(other) == type(self):
            if (self.year == other.year and self.month == other.month and
                    self.day == other.day and self.hour == other.hour and
                    self.timezone == other.timezone):
                return True
        return False

    def sameHourOfYear(self, other):
        if type(other) == type(self):
            if (self.month == other.month and self.day == other.day and
                    self.hour == other.hour and self.timezone == other.timezone):
                return True
        return False

    def hourOfYear(self):
        daysInPriorMonths = sum(self.daysPerMonth[:(self.month - 1)])
        hoursInPriorMonths = daysInPriorMonths * 24
        hoursThisMonth = (self.day - 1) * 24 + self.hour
        return hoursInPriorMonths + hoursThisMonth

    def convertToTimezone(self, tgtTimezone):
        (year, month, day, hour) = (self.year, self.month, self.day, self.hour)
        if tgtTimezone == 'EST' and self.timezone == 'UTC':
            timeChange = -4
        elif tgtTimezone == 'EST' and self.timezone == 'CST':
            timeChange = 1
        elif tgtTimezone == 'EST' and self.timezone == 'EST':
            timeChange = 0
        if timeChange != 0:
            if timeChange < 0:
                if hour > abs(timeChange):
                    hour += timeChange
                else:
                    hour = hour + 24 + timeChange
                    if day > 1:
                        day -= 1
                    else:
                        day = self.daysPerMonth[month - 1 - 1]  # shift for idx then time change
                        if month > 1:
                            month -= 1
                        else:
                            month = 12
                            year -= 1
            else:
                if hour + timeChange <= 24:
                    hour += timeChange
                else:
                    hour = hour - 24 + timeChange
                    if day < self.daysPerMonth[month - 1]:
                        day += 1
                    else:
                        day = 1
                        if month < 12:
                            month += 1
                        else:
                            month = 1
                            year += 1
        return Date(year, month, day, hour, tgtTimezone)

    def __str__(self):
        return (str(self.month) + '/' + str(self.day) + '/' + str(self.year) + ' '
                + str(self.hour) + ' ' + self.timezone)

    def __repr__(self):
        return (str(self.month) + '/' + str(self.day) + '/' + str(self.year) + ' '
                + str(self.hour) + ' ' + self.timezone)


def testDateClass():
    print('Testing date class')
    date1 = Date(2016, 5, 3, 2, 'UTC')
    date3 = Date(2016, 5, 3, 2, 'UTC')
    assert (date1 == date3)
    print(date1)
    dates = [date1, date3]
    print(dates)
    assert (date1 in dates)
    date2 = Date(2015, 5, 3, 2, 'UTC')
    date4 = Date(2016, 5, 1, 2, 'UTC')
    assert (date1 != date2)
    assert (date4 not in dates)
    assert (date1.sameHourOfYear(date2))
    assert (date1.sameHourOfYear(date3))
    assert (not date1.sameHourOfYear(date4))
    datestr = Date('2016', '5', '3', '2', 'UTC')
    assert (date1 == datestr)
    assert (dates.index(date1) == 0)
    dates.append(date4)
    assert (dates.index(date4) == 2)
    hour1 = Date(2016, 1, 1, 1, 'UTC')
    hour8760 = Date(2016, 12, 31, 24, 'UTC')
    assert (hour1.hourOfYear() == 1)
    assert (hour8760.hourOfYear() == 8760)
    assert (Date(2015, 1, 1, 1, 'UTC').hourOfYear() == 1)
    assert (Date(2015, 2, 1, 5, 'UTC').hourOfYear() == (31 * 24 + 5))
    assert (Date(2015, 2, 5, 23, 'UTC').hourOfYear() == (31 * 24 + 4 * 24 + 23))
    assert (Date(2015, 11, 27, 10, 'UTC').hourOfYear() == (8760 - 24 * 31 - 3 * 24 - 14))
    hour1est = Date(2016, 1, 1, 1, 'EST')
    assert (not hour1est == hour1)
    assert (Date(2016, 5, 3, 2, 'EST') not in dates)
    assert (date1.convertToTimezone('EST') == Date(2016, 5, 2, 22, 'EST'))
    assert (Date(2016, 5, 5, 5, 'UTC').convertToTimezone('EST') == Date(2016, 5, 5, 1, 'EST'))
    assert (Date(2016, 1, 1, 1, 'UTC').convertToTimezone('EST') == Date(2015, 12, 31, 21, 'EST'))
    assert (Date(2016, 2, 1, 1, 'UTC').convertToTimezone('EST') == Date(2016, 1, 31, 21, 'EST'))
    date5 = Date(2016, 5, 3, '00:00', 'UTC')
    date6 = Date(2016, 5, 3, '13:35', 'UTC')
    date7 = Date(2016, 5, 3, '13:55', 'UTC')
    date8 = Date(2016, 5, 3, '1:15', 'UTC')
    print(date5, date6, date7, date8)
    assert (date8 == date3)
    date9 = Date(2016, 5, 3, 2, 'CST')
    date10 = Date(2016, 5, 3, 2, 'EST')
    assert (date10.convertToTimezone('EST') == Date(2016, 5, 3, 2, 'EST'))
    assert (date10.convertToTimezone('EST') == date10)
    assert (date9.convertToTimezone('EST') == Date(2016, 5, 3, 3, 'EST'))
    assert (Date(2015, 12, 31, 24, 'CST').convertToTimezone('EST') == Date(2016, 1, 1, 1, 'EST'))
    assert (Date(2015, 1, 31, 24, 'CST').convertToTimezone('EST') == Date(2015, 2, 1, 1, 'EST'))
    date11 = Date(16, 5, 3, 2, 'EST')
    print(date11)
    assert (date11 == date10)

# testDateClass()
