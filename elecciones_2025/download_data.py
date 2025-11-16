import asyncio
from playwright.async_api import async_playwright

async def download_data():
    async with async_playwright() as p:
        browser = await p.firefox.launch(headless=False)
        ctx = await browser.new_context(accept_downloads=True)
        page = await ctx.new_page()
        await page.goto("https://elecciones.servel.cl/")

        await page.get_by_role("button", name="Diputados").click()
        await page.click("text=División Electoral Chile")
        
        for distrito in range(1, 29):
            select = page.locator("select.combo.text-filtros").nth(2)
            await select.select_option(label=f"DISTRITO {distrito}")
            
            async with page.expect_download() as download_info:
                await page.locator('a[title="Descargar en EXCEL"]:visible').click()
            
            download = await download_info.value
            await download.save_as(f"datos_2025/distrito_{distrito}.xlsx")

if __name__ == "__main__":
    asyncio.run(download_data())