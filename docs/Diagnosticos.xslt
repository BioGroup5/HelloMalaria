<?xml version="1.0" encoding="utf-8"?>
<xsl:stylesheet version="1.0" xmlns:xsl="http://www.w3.org/1999/XSL/Transform">
	<xsl:template match="/">
		<html>
			<head>
				<title>Diagnóstico de Malaria</title>
				<link rel="stylesheet" type="text/css" href="Diagnosticos.css"/>
			</head>
			<body>
				<h2>Diagnósticos de Malaria</h2>
				<div class="table-container">
				<table>
					<tr>
						<th>ID</th>
						<th>Datos personales</th>
						<th>Estado</th>
						<th>Parásito</th>
						<th>Densidad(Parasitos/μL)</th>
						<th>Episodio</th>
						<th>Infección</th>
						<th>Fecha</th>
						<th>Hora</th>
						<th>Frotis</th>
						<th>Observaciones</th>
					</tr>

					<xsl:for-each select="diagnosticos/paciente">
						<tr>
							<td><xsl:value-of select="@id"/></td>
							<td>
							  Nombre y apellido: 
							  <xsl:value-of select="datosPersonales/nombre"/>
							  <xsl:text> </xsl:text>
							  <xsl:value-of select="datosPersonales/apellido"/>
							  <br/>
							  Edad: 
							  <xsl:value-of select="datosPersonales/edad"/> años
							  <br/>
							  Género: 
							  <xsl:value-of select="datosPersonales/genero"/>
							</td>
							<td>
								<xsl:choose>
									<xsl:when test="diagnostico/estado = 'Positivo'">
										<span class="positivo">POSITIVO</span>
									</xsl:when>
									<xsl:otherwise>
										<span class="negativo">NEGATIVO</span>
									</xsl:otherwise>
								</xsl:choose>
							</td>

							<!-- Si el estado es negativo, se ocultan las demás celdas médicas -->
							<xsl:choose>
								<xsl:when test="diagnostico/estado = 'Negativo'">
									<td colspan="4" class="negativo" style="text-align:center;" >
										Sin datos clínicos relevantes
									</td>
								</xsl:when>
								<xsl:otherwise>
								  <td class="datos-positivo"><xsl:value-of select="diagnostico/parasito"/></td>
								  <td class="datos-positivo"><xsl:value-of select="diagnostico/densidadParasitaria"/></td>
								  <td class="datos-positivo"><xsl:value-of select="diagnostico/episodio"/></td>
								  <td class="datos-positivo"><xsl:value-of select="diagnostico/nivelInfeccion"/></td>
								</xsl:otherwise>
							</xsl:choose>
							<td><xsl:value-of select="fecha"/></td>
							<td><xsl:value-of select="hora"/></td>	
							<!-- Dentro de la celda de la imagen en tu XSLT -->
							<td style="text-align:center;">
							  <img src="{imagen}" alt="Imagen" style="max-width:100px; max-height:80px; cursor:pointer;" onclick="openModal(this.src)"/>
							</td>
							<td><xsl:value-of select="observaciones"/></td>
						</tr>
					</xsl:for-each>
				</table>
				</div>
			</body>
		</html>
	</xsl:template>
</xsl:stylesheet>
