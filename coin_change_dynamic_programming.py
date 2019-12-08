def application():
    change = int(input('Please, specify the amount of change: '))
    valuesUser = str(input('Please, specify the available coin denominations, separated by spaces: ')).split()
    values = [int(x) for x in valuesUser]
    def changeCoins(values, change):
        solutions = [0] * (change + 1)
        usedCoins = [0] * (change + 1)
        for coins in range(change + 1):
            minimum_solution = coins #assume that it is a change of solely 1-cent coins
            current_coin = 1 #yes, solely 1-cent coins!
            for n in [x for x in values if x <= coins]:
                if solutions[coins - n] + 1 < minimum_solution: #looking for a minimum, which is less than previous minimum
                    minimum_solution = solutions[coins - n] + 1
                    current_coin = n
            solutions[coins] = minimum_solution
            usedCoins[coins] = current_coin

        current_coins = change
        list_of_coins = []
        while current_coins > 0:
            this_coin = usedCoins[current_coins]
            list_of_coins.append(this_coin)
            current_coins = current_coins - this_coin

        print('You have requested to give a change of', change, 'coins')
        print('To do this, you need', solutions[change], 'coins:', *list_of_coins)
    changeCoins(values, change)
    
application()
