install.packages("highcharter")
install.packages("gapminder")
# Packages for creating the chart
library(highcharter)
# Other supporting packages
library(dplyr)
# Sample data for our charts
library(gapminder)
# htmlwidgets   
library(htmlwidgets)
#webshot
library(webshot)
head(gapminder)
gapminder_2007 <- gapminder:: gapminder %>% 
  filter(year == max(year))
dout_scatter <- gapminder_2007 %>% 
  hchart("scatter", hcaes(x = gdpPercap, y = lifeExp, group = country)) %>% 
  hc_title(text = "Life Expectancy vs GDP per Capita (2007)") %>% 
  hc_subtitle(text = "Data from Gapminder") %>% 
  hc_xAxis(title = list(text = "GDP per Capita")) %>% 
  hc_yAxis(title = list(text = "Life Expectancy")) %>% 
  hc_tooltip(pointFormat = "{point.country}: <b>{point.y}</b> years") %>% 
  hc_legend(enabled = FALSE)
dout_sunburst <- gapminder_2007 %>% 
  hchart("sunburst", hcaes(name = country, value = pop)) %>% 
  hc_title(text = "Population by Country (2007)") %>% 
  hc_subtitle(text = "Data from Gapminder") %>% 
  hc_tooltip(pointFormat = "<b>{point.name}</b>: {point.value}") %>% 
  hc_legend(enabled = FALSE)
dout_sunburst
dout <- data_to_hierarchical(gapminder_2007, c("continent", "country"), "pop") # explain this line of commands 
dout
hchart(dout, type = "sunburst")

pl3 <- highchart() %>%
  hc_chart(type = "sunburst") %>%
  hc_title(text = "Using Highcharter for creating a sunburst chart") %>%
  hc_subtitle(text = "#techanswers88",style = list(fontSize = "16px",fontWeight = "bold", color = "red")) %>%
  hc_add_series(name = "pop" , data = dout)%>%
  hc_plotOptions(series = list(dataLabels = list(format = "{point.name}"))) %>%
  hc_add_theme(hc_theme_economist())

#hc_theme_538

pl4 <- highchart() %>%
  hc_chart(type = "sunburst") %>%
  hc_title(text = "Using Highcharter for creating a sunburst chart") %>%
  hc_subtitle(text = "#techanswers88",style = list(fontSize = "16px",fontWeight = "bold", color = "red")) %>%
  hc_add_series(name = "pop" , data = dout)%>%
  hc_add_theme(hc_theme_538())
pl4
#hc_theme_alone()
p15 <- highchart() %>%
  hc_chart(type = "sunburst") %>%
  hc_title(text = "Using Highcharter for creating a sunburst chart") %>%
  hc_subtitle(text = "#techanswers88",style = list(fontSize = "16px",fontWeight = "bold", color = "red")) %>%
  hc_add_series(name = "pop" , data = dout)%>%
  hc_add_theme(hc_theme_alone())
p15





#Save your chart
# As highcharter charts are rendered slowly on screen so give it a delay before saving the chart.
# In our example we have used delay = 4 so that the chart gets saved after 4 seconds when it start rendering.
# If your chart is complicated and takes more time to render and your saved image is not complete then increase the delay
# eg. delay = 10
# Remember that it will save the chart which was created last. 
# If you are creating 10 charts then you can use the same command after
# creating each chart to save them one by one.
## save as a static image to be use in your word document or presentation
library(htmlwidgets)
library(webshot)

# Use these commands for saving the charts. Change the name of the plot as needed.
htmlwidgets::saveWidget(widget = pl3, file = "sunburst_chart.html", selfcontained = TRUE)
webshot::webshot("sunburst_chart.html", file = "sunburst_chart.png", delay = 4)
